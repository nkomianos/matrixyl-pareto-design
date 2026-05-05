"""
Evolutionary algorithms for peptide optimization.
Implements both single-objective (CMA-ES) and multi-objective (NSGA-II) optimization.
"""

import numpy as np
from typing import List, Tuple, Callable
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.termination import get_termination
import logging

logger = logging.getLogger(__name__)


class PeptideSequenceVariable:
    """
    Discrete variable type for peptide sequences.
    Encodes/decodes amino acids to/from integers.
    """

    AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids
    AA_TO_INT = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
    INT_TO_AA = {i: aa for i, aa in enumerate(AMINO_ACIDS)}

    @staticmethod
    def encode_sequence(sequence: str) -> np.ndarray:
        """Convert sequence string to integer array."""
        return np.array([PeptideSequenceVariable.AA_TO_INT[aa] for aa in sequence])

    @staticmethod
    def decode_sequence(encoded: np.ndarray) -> str:
        """Convert integer array back to sequence string."""
        return ''.join(PeptideSequenceVariable.INT_TO_AA[int(x)] for x in encoded)


class PeptideOptimizationProblem(Problem):
    """
    Multi-objective optimization problem for peptide design.

    Objectives:
    1. Maximize penetration score (TPSA, MW, LogP)
    2. Maximize collagen-binding affinity
    3. Minimize synthesis complexity (optional)
    """

    def __init__(
        self,
        oracle,
        constraints,
        reference_sequence: str,
        objectives: List[str] = None,
    ):
        """
        Args:
            oracle: EnsembleOracle instance
            constraints: PhysicochemicalCalculator instance
            reference_sequence: starting sequence (e.g., Matrixyl)
            objectives: ['penetration', 'binding', 'synthesis_cost'] or subset
        """
        self.oracle = oracle
        self.constraints = constraints
        self.reference = reference_sequence
        self.ref_length = len(reference_sequence)
        self.objectives = objectives or ['penetration', 'binding']

        # Variable bounds: allow up to ±2 residues
        self.min_length = max(self.ref_length - 2, 4)
        self.max_length = self.ref_length + 2

        # Problem definition: length is variable, so we fix to max and allow padding
        self.sequence_length = self.max_length

        n_vars = self.sequence_length
        n_objs = len(self.objectives)

        super().__init__(
            n_var=n_vars,
            n_obj=n_objs,
            n_constr=1,  # Constraint: length must be in valid range
            type_var=int,
            elementwise_evaluation=True
        )

    def _evaluate(self, x, out, *args, **kwargs):
        """
        Evaluate fitness for a single sequence.

        Args:
            x: integer-encoded sequence (with padding)
        """
        # Decode and handle padding
        sequence = PeptideSequenceVariable.decode_sequence(x)
        sequence = sequence.replace('A', '')  # Use 'A' as padding (can be removed)
        if not sequence or len(sequence) < self.min_length:
            # Invalid: too short
            f = [0.0] * len(self.objectives)
            g = [1.0]  # Constraint violated
            out['F'] = f
            out['G'] = g
            return

        # Ensure valid structure
        if not self.oracle.validate_structure(sequence):
            f = [0.0] * len(self.objectives)
            g = [1.0]
            out['F'] = f
            out['G'] = g
            return

        # Compute objectives
        f = []

        # Objective 1: Penetration
        if 'penetration' in self.objectives:
            penetration = self.constraints.compute_penetration_score(sequence)
            f.append(penetration)

        # Objective 2: Binding affinity
        if 'binding' in self.objectives:
            binding_scores, _ = self.oracle.score_binding_affinity([sequence])
            binding = binding_scores[0]
            f.append(binding)

        # Objective 3: Synthesis complexity (optional)
        if 'synthesis_cost' in self.objectives:
            # Simple heuristic: penalize sequences with rare residues
            rare_residues = sum(1 for aa in sequence if aa in 'CWP')
            synthesis_cost = 1.0 - (rare_residues / len(sequence))
            f.append(synthesis_cost)

        # Constraint: length in valid range
        g = [1.0 if len(sequence) < self.min_length or len(sequence) > self.max_length else 0.0]

        out['F'] = f
        out['G'] = g


class EvolutionaryOptimizer:
    """
    High-level interface for peptide optimization.
    """

    def __init__(self, oracle, constraints, reference_sequence: str):
        self.oracle = oracle
        self.constraints = constraints
        self.reference = reference_sequence
        self.history = []

    def optimize_nsga2(
        self,
        objectives: List[str] = None,
        population_size: int = 50,
        max_generations: int = 100,
    ) -> Tuple[List[str], np.ndarray]:
        """
        Multi-objective optimization using NSGA-II.

        Args:
            objectives: ['penetration', 'binding', 'synthesis_cost']
            population_size: population size
            max_generations: number of generations

        Returns:
            (pareto_sequences, objectives_matrix)
            pareto_sequences: list of Pareto-optimal sequences
            objectives_matrix: shape (n_pareto, n_objectives)
        """
        logger.info(f"Starting NSGA-II with pop={population_size}, gen={max_generations}")

        problem = PeptideOptimizationProblem(
            self.oracle,
            self.constraints,
            self.reference,
            objectives=objectives or ['penetration', 'binding']
        )

        algorithm = NSGA2(pop_size=population_size)
        termination = get_termination("n_gen", max_generations)

        res = minimize(
            problem,
            algorithm,
            termination,
            seed=42,
            verbose=True
        )

        # Extract solutions
        X = res.X  # Design variables
        F = res.F  # Objectives

        pareto_sequences = []
        for x in X:
            seq = PeptideSequenceVariable.decode_sequence(x).replace('A', '')
            if self.oracle.validate_structure(seq) and len(seq) >= 4:
                pareto_sequences.append(seq)

        self.history.append({
            'algorithm': 'NSGA2',
            'objectives': problem.objectives,
            'pareto_front': F,
            'sequences': pareto_sequences
        })

        logger.info(f"Found {len(pareto_sequences)} Pareto-optimal sequences")
        return pareto_sequences, F

    def mutate_sequence(self, sequence: str, mutation_rate: float = 0.2) -> str:
        """
        Apply random mutations to a sequence.

        Args:
            sequence: input sequence
            mutation_rate: probability of mutation per position

        Returns:
            mutated sequence
        """
        mutated = list(sequence)

        for i in range(len(mutated)):
            if np.random.random() < mutation_rate:
                # Three types of mutations
                rand = np.random.random()

                if rand < 0.6:  # Substitution (60%)
                    mutated[i] = np.random.choice(list(PeptideSequenceVariable.AMINO_ACIDS))

                elif rand < 0.8:  # Insertion (20%)
                    new_aa = np.random.choice(list(PeptideSequenceVariable.AMINO_ACIDS))
                    mutated.insert(i, new_aa)

                else:  # Deletion (20%)
                    if len(mutated) > 4:
                        mutated.pop(i)

        return ''.join(mutated)

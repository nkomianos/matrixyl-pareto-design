"""
Peptide Optimization Framework
Core modules for evolutionary optimization of cosmetic peptides.
"""

from .candidates import CandidateEvaluation, CandidateEvaluator
from .analysis import CandidateRow
from .constraints import PenetrationScore, PhysicochemicalCalculator
from .chemistry import MolecularDescriptors
from .experimental_design import CompactnessRow, SynthesisCandidate
from .function_scores import FunctionalPreservationScore, FunctionalPreservationScorer
from .pareto_search import (
    ParetoCandidate,
    ParetoGenerationSummary,
    ParetoObjectives,
    ParetoSearchConfig,
    ParetoSearchOptimizer,
    ParetoSearchResult,
)
from .search_space import SearchSpaceRules, SearchSpaceValidation
from .structure import ConformerEnsembleSummary, ConformerMetrics
from .tournament_search import (
    GenerationSummary,
    TournamentSearchConfig,
    TournamentSearchOptimizer,
    TournamentSearchResult,
)

__version__ = '0.1.0'
__all__ = [
    'CandidateEvaluation',
    'CandidateEvaluator',
    'CandidateRow',
    'CompactnessRow',
    'FunctionalPreservationScore',
    'FunctionalPreservationScorer',
    'MolecularDescriptors',
    'ParetoCandidate',
    'ParetoGenerationSummary',
    'ParetoObjectives',
    'ParetoSearchConfig',
    'ParetoSearchOptimizer',
    'ParetoSearchResult',
    'PenetrationScore',
    'PhysicochemicalCalculator',
    'SearchSpaceRules',
    'SearchSpaceValidation',
    'SynthesisCandidate',
    'ConformerEnsembleSummary',
    'ConformerMetrics',
    'GenerationSummary',
    'TournamentSearchConfig',
    'TournamentSearchOptimizer',
    'TournamentSearchResult',
]

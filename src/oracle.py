"""
Protein Language Model Oracle
Wraps multiple pretrained PLMs (ESM-2, ProtBERT) for sequence scoring.
Used as fitness function oracle in evolutionary algorithm.
"""

import torch
import numpy as np
from typing import Tuple, List
from transformers import AutoTokenizer, AutoModel
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class EnsembleOracle:
    """
    Ensemble of protein language models for binding affinity prediction.

    Attributes:
        models: dict of loaded model instances
        tokenizers: dict of tokenizers
        device: torch device (cuda/cpu)
        uncertainty_estimation: if True, returns (score, std_dev) tuple
    """

    def __init__(
        self,
        models: List[str] = None,
        device: str = None,
        batch_size: int = 32,
        uncertainty_estimation: bool = True
    ):
        """
        Initialize ensemble oracle with pretrained models.

        Args:
            models: List of model names. Default: ['esm2_t33', 'protbert']
            device: 'cuda' or 'cpu'. Auto-detect if None.
            batch_size: Inference batch size
            uncertainty_estimation: Compute std dev across ensemble
        """
        if models is None:
            models = ['facebook/esm2_t33_650M_UR50D', 'rostlab/prot_bert']

        self.device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
        self.batch_size = batch_size
        self.uncertainty_estimation = uncertainty_estimation
        self.models = {}
        self.tokenizers = {}

        logger.info(f"Initializing oracle on {self.device}")
        for model_name in models:
            logger.info(f"Loading {model_name}...")
            self.tokenizers[model_name] = AutoTokenizer.from_pretrained(model_name)
            self.models[model_name] = AutoModel.from_pretrained(
                model_name,
                output_hidden_states=True
            ).to(self.device).eval()

    def _get_embeddings(self, sequences: List[str], model_name: str) -> np.ndarray:
        """
        Get sequence embeddings from a single model.
        Returns mean pooling over sequence length.
        """
        tokenizer = self.tokenizers[model_name]
        model = self.models[model_name]

        embeddings = []
        with torch.no_grad():
            for i in range(0, len(sequences), self.batch_size):
                batch = sequences[i:i + self.batch_size]

                # Tokenize
                tokens = tokenizer(
                    batch,
                    return_tensors='pt',
                    padding=True,
                    truncation=True,
                    max_length=1024
                ).to(self.device)

                # Inference
                outputs = model(**tokens)

                # Mean pooling over sequence dimension
                # Shape: (batch_size, hidden_size)
                last_hidden = outputs.last_hidden_state
                mask = tokens['attention_mask'].unsqueeze(-1).expand(last_hidden.size()).float()
                sum_hidden = (last_hidden * mask).sum(dim=1)
                seq_length = mask.sum(dim=1)
                mean_pooled = sum_hidden / seq_length

                embeddings.append(mean_pooled.cpu().numpy())

        return np.vstack(embeddings)

    def score_binding_affinity(
        self,
        sequences: List[str] or str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Score collagen-binding affinity via ensemble embedding similarity.

        Args:
            sequences: Single sequence or list of sequences

        Returns:
            (scores, uncertainties): arrays of shape (n_seqs,)
            scores: ensemble mean affinity score [0, 1]
            uncertainties: std dev across ensemble models
        """
        if isinstance(sequences, str):
            sequences = [sequences]

        # Reference collagen-binding motif (from literature on Matrixyl)
        # This is a simplified proxy; in practice, could use in vitro binding data
        reference_seq = "GPKGDP"  # Collagen-targeting tripeptide core

        all_scores = []
        for model_name in self.models.keys():
            # Get embeddings for all sequences + reference
            test_embeddings = self._get_embeddings(sequences, model_name)
            ref_embedding = self._get_embeddings([reference_seq], model_name)[0]

            # Cosine similarity
            from sklearn.metrics.pairwise import cosine_similarity
            similarities = cosine_similarity(test_embeddings, [ref_embedding]).flatten()

            # Normalize to [0, 1]
            scores = (similarities + 1) / 2  # cosine similarity is [-1, 1]
            all_scores.append(scores)

        all_scores = np.array(all_scores)  # shape: (n_models, n_seqs)

        ensemble_mean = all_scores.mean(axis=0)
        ensemble_std = all_scores.std(axis=0)

        return ensemble_mean, ensemble_std

    def validate_structure(self, sequence: str) -> bool:
        """
        Check if sequence is likely to fold properly (not intrinsically disordered).
        Simple heuristic: check amino acid composition.

        Returns:
            True if sequence passes basic validity checks
        """
        if len(sequence) < 4 or len(sequence) > 100:
            return False

        # Reject sequences with >30% disorder-promoting residues (P, E, S alone)
        disorder_residues = sum(1 for aa in sequence if aa in 'PES')
        if disorder_residues / len(sequence) > 0.3:
            return False

        # Reject sequences with rare amino acids
        valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(aa in valid_aas for aa in sequence):
            return False

        return True

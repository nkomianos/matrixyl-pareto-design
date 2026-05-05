"""
Protein language model oracle for Matrixyl-family candidate scoring.

The module is import-safe on machines without GPU dependencies. Heavy
dependencies are loaded only when an EnsembleOracle instance is created.
"""

from dataclasses import asdict, dataclass
import hashlib
import logging
from pathlib import Path
import re
from typing import Iterable, Optional

import numpy as np


logger = logging.getLogger(__name__)

DEFAULT_REFERENCE_SEQUENCE = "KTTKS"
DEFAULT_MODEL_NAMES = (
    "facebook/esm2_t33_650M_UR50D",
    "rostlab/prot_bert",
)
CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
PLM_ALLOWED_AMINO_ACIDS = CANONICAL_AMINO_ACIDS | frozenset("UZOBX")
UNCOMMON_AA_PATTERN = re.compile(r"[UZOB]")


class PLMDependencyError(ImportError):
    """Raised when PLM scoring is requested without torch/transformers."""


@dataclass(frozen=True)
class ProteinModelSpec:
    model_name: str
    input_format: str = "raw"
    max_length: int = 1024
    do_lower_case: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


def normalize_plm_sequence(sequence: str) -> str:
    normalized = "".join(sequence.strip().upper().split())
    if not normalized:
        raise ValueError("sequence cannot be empty")
    invalid = sorted(set(normalized) - PLM_ALLOWED_AMINO_ACIDS)
    if invalid:
        raise ValueError(f"invalid protein sequence characters: {''.join(invalid)}")
    return normalized


def model_spec_from_name(model_name: str) -> ProteinModelSpec:
    lower = model_name.lower()
    if "prot_bert" in lower or "protbert" in lower:
        return ProteinModelSpec(model_name=model_name, input_format="spaced")
    return ProteinModelSpec(model_name=model_name, input_format="raw")


def model_specs_from_names(model_names: Iterable[str]) -> tuple[ProteinModelSpec, ...]:
    return tuple(model_spec_from_name(model_name) for model_name in model_names)


def format_sequence_for_model(sequence: str, spec: ProteinModelSpec) -> str:
    normalized = normalize_plm_sequence(sequence)
    if spec.input_format == "spaced":
        protbert_safe = UNCOMMON_AA_PATTERN.sub("X", normalized)
        return " ".join(protbert_safe)
    if spec.input_format == "raw":
        return normalized
    raise ValueError(f"Unsupported model input format: {spec.input_format}")


def cosine_similarity_matrix(embeddings: np.ndarray, reference_embedding: np.ndarray) -> np.ndarray:
    if embeddings.ndim != 2:
        raise ValueError("embeddings must be a 2D array")
    reference = np.asarray(reference_embedding).reshape(1, -1)
    numerator = embeddings @ reference.T
    denominator = np.linalg.norm(embeddings, axis=1, keepdims=True) * np.linalg.norm(reference)
    denominator = np.where(denominator == 0, 1.0, denominator)
    return (numerator / denominator).reshape(-1)


class EmbeddingCache:
    """Small disk cache for per-model, per-sequence embedding vectors."""

    def __init__(self, cache_dir: str | Path):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _safe_model_name(model_name: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", model_name)

    def path_for(self, model_name: str, sequence: str) -> Path:
        digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]
        safe_model = self._safe_model_name(model_name)
        return self.cache_dir / safe_model / f"{sequence}_{digest}.npy"

    def get(self, model_name: str, sequence: str) -> Optional[np.ndarray]:
        path = self.path_for(model_name, sequence)
        if not path.exists():
            return None
        return np.load(path)

    def set(self, model_name: str, sequence: str, embedding: np.ndarray) -> None:
        path = self.path_for(model_name, sequence)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path, np.asarray(embedding, dtype=np.float32))


def _load_plm_dependencies():
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:
        raise PLMDependencyError(
            "PLM scoring requires torch and transformers. Install requirements.txt "
            "on the GPU host before running the oracle."
        ) from exc
    return torch, AutoModel, AutoTokenizer


class EnsembleOracle:
    """
    Ensemble of protein language models for Matrixyl-family similarity scoring.

    Scores are embedding cosine similarity to a configurable reference sequence,
    normalized from [-1, 1] to [0, 1]. This is a computational proxy for
    functional preservation, not a direct biological activity assay.
    """

    def __init__(
        self,
        models: Optional[list[str]] = None,
        *,
        reference_sequence: str = DEFAULT_REFERENCE_SEQUENCE,
        device: Optional[str] = None,
        batch_size: int = 16,
        cache_dir: Optional[str | Path] = None,
        local_files_only: bool = False,
    ):
        if batch_size < 1:
            raise ValueError("batch_size must be at least 1")

        self.torch, AutoModel, AutoTokenizer = _load_plm_dependencies()
        self.reference_sequence = normalize_plm_sequence(reference_sequence)
        self.model_specs = model_specs_from_names(models or DEFAULT_MODEL_NAMES)
        self.device = device or ("cuda" if self.torch.cuda.is_available() else "cpu")
        self.batch_size = batch_size
        self.cache = EmbeddingCache(cache_dir) if cache_dir else None
        self.models = {}
        self.tokenizers = {}

        logger.info("Initializing PLM oracle on %s", self.device)
        for spec in self.model_specs:
            logger.info("Loading %s", spec.model_name)
            tokenizer_kwargs = {"local_files_only": local_files_only}
            if spec.input_format == "spaced":
                tokenizer_kwargs["do_lower_case"] = spec.do_lower_case
            self.tokenizers[spec.model_name] = AutoTokenizer.from_pretrained(
                spec.model_name,
                **tokenizer_kwargs,
            )
            self.models[spec.model_name] = AutoModel.from_pretrained(
                spec.model_name,
                output_hidden_states=True,
                local_files_only=local_files_only,
            ).to(self.device).eval()

    def _embed_batch(self, sequences: list[str], spec: ProteinModelSpec) -> np.ndarray:
        tokenizer = self.tokenizers[spec.model_name]
        model = self.models[spec.model_name]
        formatted = [format_sequence_for_model(sequence, spec) for sequence in sequences]

        tokens = tokenizer(
            formatted,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=spec.max_length,
            return_special_tokens_mask=True,
        )
        special_tokens_mask = tokens.pop("special_tokens_mask", None)
        tokens = {key: value.to(self.device) for key, value in tokens.items()}
        if special_tokens_mask is not None:
            special_tokens_mask = special_tokens_mask.to(self.device)
        with self.torch.no_grad():
            outputs = model(**tokens)
            last_hidden = outputs.last_hidden_state
            token_mask = tokens["attention_mask"]
            if special_tokens_mask is not None:
                token_mask = token_mask * (1 - special_tokens_mask)
            mask = token_mask.unsqueeze(-1).expand(last_hidden.size()).float()
            summed = (last_hidden * mask).sum(dim=1)
            lengths = mask.sum(dim=1).clamp(min=1.0)
            pooled = summed / lengths
        return pooled.detach().cpu().numpy()

    def _get_embeddings(self, sequences: list[str], spec: ProteinModelSpec) -> np.ndarray:
        normalized = [normalize_plm_sequence(sequence) for sequence in sequences]
        embeddings: list[Optional[np.ndarray]] = []
        missing_sequences: list[str] = []
        missing_indices: list[int] = []

        for index, sequence in enumerate(normalized):
            cached = self.cache.get(spec.model_name, sequence) if self.cache else None
            embeddings.append(cached)
            if cached is None:
                missing_sequences.append(sequence)
                missing_indices.append(index)

        for start in range(0, len(missing_sequences), self.batch_size):
            batch = missing_sequences[start:start + self.batch_size]
            batch_embeddings = self._embed_batch(batch, spec)
            if len(batch_embeddings) != len(batch):
                raise RuntimeError(
                    f"Model {spec.model_name} returned {len(batch_embeddings)} embeddings "
                    f"for a batch of {len(batch)} sequences"
                )
            for sequence, embedding, original_index in zip(
                batch,
                batch_embeddings,
                missing_indices[start:start + self.batch_size],
            ):
                if self.cache:
                    self.cache.set(spec.model_name, sequence, embedding)
                embeddings[original_index] = embedding

        if any(embedding is None for embedding in embeddings):
            raise RuntimeError(f"Missing embeddings after scoring with {spec.model_name}")
        return np.vstack(embeddings)

    def score_binding_affinity(
        self,
        sequences: list[str] | str,
        *,
        reference_sequence: Optional[str] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(sequences, str):
            sequences = [sequences]
        normalized = [normalize_plm_sequence(sequence) for sequence in sequences]
        reference = normalize_plm_sequence(reference_sequence or self.reference_sequence)

        model_scores = []
        for spec in self.model_specs:
            candidate_embeddings = self._get_embeddings(normalized, spec)
            reference_embedding = self._get_embeddings([reference], spec)[0]
            similarities = cosine_similarity_matrix(candidate_embeddings, reference_embedding)
            model_scores.append((similarities + 1.0) / 2.0)

        scores = np.vstack(model_scores)
        return scores.mean(axis=0), scores.std(axis=0)

    def validate_structure(self, sequence: str) -> bool:
        normalized = normalize_plm_sequence(sequence)
        if len(normalized) < 4 or len(normalized) > 100:
            return False
        disorder_residues = sum(1 for aa in normalized if aa in "PES")
        return disorder_residues / len(normalized) <= 0.6

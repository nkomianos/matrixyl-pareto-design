import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from src.oracle import (
    EmbeddingCache,
    EnsembleOracle,
    ProteinModelSpec,
    cosine_similarity_matrix,
    format_sequence_for_model,
    model_spec_from_name,
    normalize_plm_sequence,
)


class OraclePrepTests(unittest.TestCase):
    def test_model_specs_capture_protbert_formatting(self):
        protbert = model_spec_from_name("Rostlab/prot_bert")
        esm = model_spec_from_name("facebook/esm2_t33_650M_UR50D")

        self.assertEqual(protbert.input_format, "spaced")
        self.assertEqual(esm.input_format, "raw")
        self.assertEqual(format_sequence_for_model("KTTKS", protbert), "K T T K S")
        self.assertEqual(format_sequence_for_model("KTTKS", esm), "KTTKS")

    def test_protbert_format_replaces_uncommon_amino_acids(self):
        spec = ProteinModelSpec("Rostlab/prot_bert", input_format="spaced")

        self.assertEqual(normalize_plm_sequence(" a u z o b x "), "AUZOBX")
        self.assertEqual(format_sequence_for_model("AUZOBX", spec), "A X X X X X")

    def test_cosine_similarity_matrix(self):
        embeddings = np.array([[1.0, 0.0], [0.0, 1.0]])
        reference = np.array([1.0, 0.0])

        similarities = cosine_similarity_matrix(embeddings, reference)

        self.assertTrue(np.allclose(similarities, np.array([1.0, 0.0])))

    def test_embedding_cache_round_trip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = EmbeddingCache(tmpdir)
            embedding = np.array([1.0, 2.0, 3.0], dtype=np.float32)

            self.assertIsNone(cache.get("model/name", "KTTKS"))
            cache.set("model/name", "KTTKS", embedding)

            cached = cache.get("model/name", "KTTKS")
            self.assertTrue(np.allclose(cached, embedding))

    def test_score_binding_affinity_can_be_exercised_without_loading_models(self):
        oracle = object.__new__(EnsembleOracle)
        oracle.reference_sequence = "KTTKS"
        oracle.model_specs = (ProteinModelSpec("dummy"),)
        oracle.batch_size = 8
        oracle.cache = None

        vectors = {
            "KTTKS": np.array([1.0, 0.0]),
            "PTTPS": np.array([0.0, 1.0]),
        }

        def embed_batch(sequences, spec):
            return np.vstack([vectors[sequence] for sequence in sequences])

        oracle._embed_batch = embed_batch
        scores, uncertainties = oracle.score_binding_affinity(["KTTKS", "PTTPS"])

        self.assertTrue(np.allclose(scores, np.array([1.0, 0.5])))
        self.assertTrue(np.allclose(uncertainties, np.array([0.0, 0.0])))

    def test_embedding_pooling_excludes_special_tokens(self):
        try:
            import torch
        except ImportError:
            self.skipTest("torch is not installed")

        class FakeTokenizer:
            def __call__(self, sequences, **kwargs):
                return {
                    "input_ids": torch.tensor([[1, 2, 3, 0]]),
                    "attention_mask": torch.tensor([[1, 1, 1, 0]]),
                    "special_tokens_mask": torch.tensor([[1, 0, 1, 0]]),
                }

        class FakeModel:
            def __call__(self, **kwargs):
                self.kwargs = kwargs
                self.last_hidden = torch.tensor(
                    [[[100.0, 100.0], [2.0, 4.0], [50.0, 50.0], [0.0, 0.0]]]
                )
                return SimpleNamespace(last_hidden_state=self.last_hidden)

        oracle = object.__new__(EnsembleOracle)
        oracle.device = "cpu"
        oracle.torch = torch
        spec = ProteinModelSpec("dummy")
        fake_model = FakeModel()
        oracle.tokenizers = {"dummy": FakeTokenizer()}
        oracle.models = {"dummy": fake_model}

        embedding = oracle._embed_batch(["KTTKS"], spec)

        self.assertNotIn("special_tokens_mask", fake_model.kwargs)
        self.assertTrue(np.allclose(embedding, np.array([[2.0, 4.0]])))

    def test_embedding_batch_size_mismatch_is_rejected(self):
        oracle = object.__new__(EnsembleOracle)
        oracle.batch_size = 8
        oracle.cache = None
        oracle._embed_batch = lambda sequences, spec: np.zeros((1, 2))

        with self.assertRaises(RuntimeError):
            oracle._get_embeddings(["KTTKS", "PTTPS"], ProteinModelSpec("dummy"))


if __name__ == "__main__":
    unittest.main()

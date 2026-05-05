import unittest


class LightweightImportTests(unittest.TestCase):
    def test_package_import_does_not_require_optimizer_dependencies(self):
        import src
        from src import PhysicochemicalCalculator
        from src.constraints import PhysicochemicalCalculator as DirectCalculator

        self.assertIs(PhysicochemicalCalculator, DirectCalculator)
        self.assertEqual(src.__version__, "0.1.0")

    def test_matrixyl_core_descriptor_shape(self):
        from src.constraints import PhysicochemicalCalculator

        calculator = PhysicochemicalCalculator()
        properties = calculator.compute_amino_acid_properties("KTTKS")

        self.assertEqual(properties["length"], 5)
        self.assertEqual(properties["charge"], 2)
        self.assertGreater(properties["mw"], 0)
        self.assertGreaterEqual(calculator.compute_penetration_score("KTTKS"), 0.0)
        self.assertLessEqual(calculator.compute_penetration_score("KTTKS"), 1.0)


if __name__ == "__main__":
    unittest.main()

import unittest
import pandas
from src.lib import ThinFilmDatabase


class TestLib(unittest.TestCase):

    def test_database(self):
        df = ThinFilmDatabase.load_df()
        self.assertIsInstance(df, pandas.DataFrame)
        rows, cols = df.shape
        self.assertEqual(rows,  8)
        self.assertEqual(cols,  8)

    def test_read_conductivity(self):
        self.assertIsInstance(
            ThinFilmDatabase.read_conductivity('TL10_10'), pandas.Series
        )

    def test_getters(self):
        self.assertIsInstance(
            ThinFilmDatabase.plt_params('TL10_10'), dict
        )


if __name__ == '__main__':
    unittest.main()

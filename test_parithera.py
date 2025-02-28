import unittest
from scanpy_scripts.parithera_umap import UMAP
from scanpy_scripts.parithera_tsne import TSNE
from scanpy_scripts.parithera_marker_genes import Markers
from scanpy_scripts.parithera_leiden import Leiden
from scanpy_scripts.parithera_cluster import Clusters
from scanpy_scripts.link_to_project import LinkToProject

class TestParithera(unittest.TestCase):
    
    output_path = "./python"

    def test_link_to_project(self):
        link_to_project = LinkToProject(TestParithera.output_path)
        import os
        file_size_mb = os.path.getsize(link_to_project.output_path.replace("python", "out.h5ad")) / (1024 * 1024)
        self.assertGreater(file_size_mb, 40, "The size of the file out.h5ad is not greater than 40MB")

    def test_umap(self):
        umap = UMAP(TestParithera.output_path)
        self.assertEqual(umap.run()["type"], 'umap', 'Err')

    def test_tsne(self):
        tsne = TSNE(TestParithera.output_path)
        self.assertEqual(tsne.run()["type"], 'tsne', 'Err')

    def test_leiden(self):
        leiden = Leiden(TestParithera.output_path)
        self.assertEqual(leiden.run()["type"], 'leiden', 'Err')

    def test_cluster(self):
        clusters = Clusters(TestParithera.output_path)
        self.assertEqual(clusters.run()["type"], 'cluster', 'Err')

    def test_marker_genes(self):
        marker = Markers(TestParithera.output_path)
        expected_genes = ["IFITM3", "A2M", "IGFBP7", "HSPG2", "IFI27"]
        output = marker.run()
        self.assertEqual(output["cluster_genes"], expected_genes, f"Expected genes {expected_genes}, but got {marker.run()['cluster_genes']}")

if __name__ == '__main__':
    unittest.main()
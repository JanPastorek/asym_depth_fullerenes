import argparse
import csv
import math
import subprocess
import tempfile
from pathlib import Path

import networkx as nx
import numpy as np


def _extract_faces_from_embedding(embedding):
    """Helper to traverse the planar embedding and return cycles of faces."""
    faces = []
    visited_edges = set()
    for u, v in embedding.edges():
        if (u, v) in visited_edges:
            continue
        try:
            face = embedding.traverse_face(u, v, visited_edges)
        except (KeyError, nx.NetworkXException):
            continue
        faces.append(face)
    return faces

def hex_neighbor_index(g6_string):
    """
    Calculates the Hexagon Neighbor Index (H) for fullerenes.
    Formula: H = sum(k^2 * h_k)
    Weights hexagons with the largest number of hexagonal neighbors most heavily.
    """
    g6_string = g6_string.strip()
    if not g6_string:
        return math.nan

    # 1. Create graph from g6
    graph = nx.from_graph6_bytes(g6_string.encode())
    
    # 2. Get Planar Embedding
    is_planar, embedding = nx.check_planarity(graph)
    if not is_planar:
        return math.nan

    # 3. Robust Face Extraction using embedding traversal
    faces = _extract_faces_from_embedding(embedding)

    # 4. Identify Hexagons (Length 6)
    # Note: Fullerenes also have 12 pentagons (Length 5)
    is_hex = [len(face) == 6 for face in faces]
    hex_face_indices = [i for i, val in enumerate(is_hex) if val]

    if not hex_face_indices:
        return 0.0

    # 5. Build Edge-to-Face Map
    # An edge in a fullerene is shared by exactly 2 faces.
    edge_to_faces = {}
    for face_idx, face in enumerate(faces):
        for i in range(len(face)):
            u, v = face[i], face[(i + 1) % len(face)]
            edge = tuple(sorted((u, v)))
            edge_to_faces.setdefault(edge, []).append(face_idx)

    # 6. Calculate k (number of hexagonal neighbors) for each hexagon
    h_index = 0
    for f_idx in hex_face_indices:
        k = 0
        face = faces[f_idx]
        for i in range(len(face)):
            u, v = face[i], face[(i + 1) % len(face)]
            edge = tuple(sorted((u, v)))
            
            # Find the other face sharing this edge
            for neighbor_face_idx in edge_to_faces.get(edge, []):
                if neighbor_face_idx != f_idx:
                    if is_hex[neighbor_face_idx]:
                        k += 1
        
        # Apply the formula H = sum(k^2 * h_k)
        # Summing k^2 for every individual hexagon is mathematically equivalent
        h_index += (k ** 2)

    return h_index


def homo_lumo_gap(g6_string):
    """
    Compute HOMO, LUMO, and the HOMO-LUMO gap from the adjacency spectrum.
    """
    g6_string = g6_string.strip()
    if not g6_string:
        return math.nan, math.nan, math.nan

    graph = nx.from_graph6_bytes(g6_string.encode())
    n = graph.number_of_nodes()
    adj = nx.to_numpy_array(graph)

    eigenvalues = np.linalg.eigvalsh(adj)
    eigenvalues.sort()

    homo_idx = (n // 2) - 1
    lumo_idx = n // 2

    homo = float(eigenvalues[homo_idx])
    lumo = float(eigenvalues[lumo_idx])
    gap = lumo - homo

    return homo, lumo, gap


def wiener_index(g6_string):
    """Compute the Wiener index: sum of shortest-path distances over all pairs."""
    g6_string = g6_string.strip()
    if not g6_string:
        return math.nan

    graph = nx.from_graph6_bytes(g6_string.encode())
    return float(nx.wiener_index(graph))


def run_asymmetric_depth_batch(julia_bin, project_dir, input_path, output_path):
    script_path = Path(__file__).resolve().parent / "asymmetric_depth_batch.jl"
    cmd = [
        julia_bin,
        f"--project={project_dir}",
        str(script_path),
        str(input_path),
        str(output_path),
    ]
    subprocess.run(cmd, check=True)


def read_g6_lines(path):
    lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                lines.append(stripped)
    return lines


def read_depths(path):
    depths = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            depths[row["g6"]] = float(row["asymmetric_depth"])
    return depths


def write_results(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "g6",
                "hex_neighbor_index",
                "homo",
                "lumo",
                "gap",
                "wiener_index",
                "asymmetric_depth",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--julia", default=None)
    parser.add_argument("--project", default=None)
    args = parser.parse_args()

    julia_bin = args.julia or "julia"
    project_dir = args.project
    if project_dir is None:
        project_dir = str(Path(__file__).resolve().parents[1])

    g6_lines = read_g6_lines(args.input)
    if not g6_lines:
        write_results(args.output, [])
        return

    rows = []
    for g6_line in g6_lines:
        hni = hex_neighbor_index(g6_line)
        homo, lumo, gap = homo_lumo_gap(g6_line)
        wiener = wiener_index(g6_line)
        rows.append({
            "g6": g6_line,
            "hex_neighbor_index": hni,
            "homo": homo,
            "lumo": lumo,
            "gap": gap,
            "wiener_index": wiener,
        })

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        input_path = tmpdir / "chunk.g6"
        output_path = tmpdir / "depths.csv"
        input_path.write_text("\n".join(g6_lines), encoding="utf-8")

        run_asymmetric_depth_batch(
            julia_bin, project_dir, input_path, output_path)
        depths = read_depths(output_path)

    for row in rows:
        row["asymmetric_depth"] = depths.get(row["g6"], math.nan)

    write_results(args.output, rows)


if __name__ == "__main__":
    main()

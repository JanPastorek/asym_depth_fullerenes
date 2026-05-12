import argparse
from pathlib import Path


def split_g6(input_path, out_dir, lines_per_chunk):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chunk_index = 0
    current_lines = []

    with open(input_path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            current_lines.append(stripped)
            if len(current_lines) >= lines_per_chunk:
                chunk_path = out_dir / f"chunk_{chunk_index:05d}.g6"
                chunk_path.write_text(
                    "\n".join(current_lines), encoding="utf-8")
                chunk_index += 1
                current_lines = []

    if current_lines:
        chunk_path = out_dir / f"chunk_{chunk_index:05d}.g6"
        chunk_path.write_text("\n".join(current_lines), encoding="utf-8")
        chunk_index += 1

    return chunk_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--lines-per-chunk", type=int, default=1000)
    args = parser.parse_args()

    count = split_g6(args.input, args.out_dir, args.lines_per_chunk)
    print(count)


if __name__ == "__main__":
    main()

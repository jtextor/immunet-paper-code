import json
import pandas as pd
from config import data_folder, cache_folder

if not data_folder.is_dir():
    data_folder.mkdir()

data = []
datasets = [f.name for f in cache_folder.iterdir() if f.is_dir()]

for dataset in datasets:

    dataset_folder = cache_folder / dataset
    slide_tiles_files = [x for x in dataset_folder.glob("*_tiles.json")]

    for slide_tiles_file in slide_tiles_files:
        slide = slide_tiles_file.stem.replace("_tiles", "")
        print(dataset, slide)

        with open(slide_tiles_file) as f:
            tiles = json.load(f)
        
        x_min = min([i["x"] for i in tiles])*2
        x_max = max([i["x"] for i in tiles])*2 + tiles[0]["width"]*2

        y_min = min([i["y"] for i in tiles])*2
        y_max = max([i["y"] for i in tiles])*2 + tiles[0]["height"]*2

        total_width = int(x_max - x_min + 0.5)
        total_height = int(y_max - y_min + 0.5)

        data.append({
            "dataset": dataset, "slide": slide,
            "x_min": x_min, "x_max": x_max,
            "y_min": y_min, "y_max": y_max
        })
        
df = pd.DataFrame(data)
df.to_csv(data_folder / "dims.csv", index=False)

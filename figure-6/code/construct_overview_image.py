import pandas as pd 
from PIL import Image
import numpy as np
from config import data_folder, images_folder

if not data_folder.is_dir():
    data_folder.mkdir()

if not images_folder.is_dir():
    images_folder.mkdir()

df = pd.read_csv(data_folder / "low_high.csv")

images = {}

max_width = 0
max_height = 0

for index, row in df.iterrows():
    dataset = row['dataset']
    slide = row['slide']
    filename = images_folder / f"{dataset}.{slide}.segmentation.png"
    print(filename)
    im = Image.open(filename)
    im = im.resize((im.width//4, im.height//4))
    f = f"{dataset}.{slide}"
    images[f] = im
    if im.width > max_width:
        max_width = im.width
    if im.height > max_height:
        max_height = im.height


print(max_width, max_height)


for low in [False, True]:
    overview = np.ones((
        max_height, 
        max_width*4,
        3), dtype=np.uint8)*255
    local_df = df[df.low==low].sort_values(by=["tissue"])
    idx = 0
    for index, row in local_df.iterrows():
        dataset = row['dataset']
        slide = row['slide']
        f = f"{dataset}.{slide}"
        im = images[f]
        x_offset = (max_width-im.width)//2
        y_offset = (max_height-im.height)//2
        arr = np.array(im)
        overview[
                y_offset:y_offset + im.height,
                x_offset + idx * max_width:x_offset+ idx * max_width + im.width
                ] = arr
        print(arr.shape)
        idx += 1
    overview_im = Image.fromarray(overview)
    if low:
        overview_im.save(images_folder / "low_density.png")
    else:
        overview_im.save(images_folder / "high_density.png")

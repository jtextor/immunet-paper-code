import numpy as np
import json
from PIL import Image
from skimage import draw
import cv2
import gzip
from utils import *

COLORS = [[255, 0, 0], [0, 255, 0], [255, 0, 255], [255, 255, 0], [0, 255, 255]]
CHANNEL_NAMES = ["cd3", "foxp3", "cd20", "cd45ro", "cd8"]
HALF_SIDE = 37

class Patch:
    def __init__(self, x_lim, y_lim):
        self.x_lim = x_lim
        self.y_lim = y_lim


def draw_ann_locations(image, annotations):

    for ann in annotations:
        iy, ix = int(round(ann["y"])), int(round(ann["x"]))
        rr, cc = draw.disk((iy, ix), 1.4, shape=image.shape)
        image[rr, cc, :] = 255

    return image


def mask_unknown(image, mask):
    mask = mask.astype(np.uint8) * 255
    mask_rgb = np.dstack((mask, mask, mask))
    masked = cv2.addWeighted(image, 1, mask_rgb, 0.5, 0)

    return masked


def extract_patches(anns, offset=HALF_SIDE):
    patches = []

    for ann in anns:
        x_lim = [ann["x"] - offset, ann["x"] + offset]
        y_lim = [ann["y"] - offset, ann["y"] + offset]
        patches.append(Patch(x_lim, y_lim))

    return patches


def save_image_patch(image, patch_xlim, patch_ylim, out_file):
    patch = image[patch_ylim[0]:patch_ylim[1], patch_xlim[0]:patch_xlim[1]]
    patch_img = Image.fromarray(patch)
    patch_img.save(out_file)


if __name__ == '__main__':

    IMAGE_FOLDER.mkdir(exist_ok=True)

    with gzip.open(ANNOTATIONS_FILE) as f:
        annotations = json.loads(f.read())

    # Annotations to center patches at
    ann1 = [ann for ann in annotations if ann["_id"] == "62fd26193de074354d0964a7"][0]
    ann2 = [ann for ann in annotations if ann["_id"] == "62fd23fd3de074354d096481"][0]

    patches = extract_patches([ann1, ann2])

    im_frame = Image.open(COMPOSITE_FILE)
    np_frame = np.array(im_frame.convert('RGB'))
    # Plot locations of annotations
    np_frame = draw_ann_locations(np_frame, annotations)

    # Generate label maps
    labels = extract_labels(np_frame, annotations)

    # Add a transparent overlay to separate regions with and without known annotation status
    mask = labels[..., 0] == -1
    masked = mask_unknown(np_frame, mask)

    # Save masked image patches
    for ind, patch in enumerate(patches):
        save_image_patch(masked, patch.x_lim, patch.y_lim, IMAGE_FOLDER / "p{}.png".format(ind + 1))

    # Save proximity maps
    proximity = (labels[..., 0] + 2) * 255 / 7
    proximity = np.round(proximity).astype(np.uint8)

    for ind, patch in enumerate(patches):
        save_image_patch(proximity, patch.x_lim, patch.y_lim, IMAGE_FOLDER / "p{}_proximity.png".format(ind + 1))

    # Save phenotype maps
    for i in range(5):
        channel = labels[..., i + 1]
        color = COLORS[i]
        chanel_name = CHANNEL_NAMES[i]

        non_zero_inds = channel > 0
        color_array = np.repeat(np.array([color]), channel[non_zero_inds].shape[0], axis=0)
        weights = np.array([channel[non_zero_inds]]).T

        channel_rbg = np.zeros(channel.shape + (3,))
        channel_rbg[non_zero_inds] = np.round(color_array * weights)
        channel_rbg = channel_rbg.astype(np.uint8).clip(0, 255)

        for ind, patch in enumerate(patches):
            save_image_patch(channel_rbg, patch.x_lim, patch.y_lim, IMAGE_FOLDER / "p{}_{}.png".format(ind + 1, chanel_name))

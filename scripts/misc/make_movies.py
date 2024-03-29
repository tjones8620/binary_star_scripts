import os
import moviepy.video.io.ImageSequenceClip
import re
import numpy as np

def make_movies(image_folder, movie_folder, movie_name):
    image_folder
    fps=10
    dir_list = os.listdir(image_folder)

    dir_list.sort(key=lambda test_string : list(
        map(int, re.findall(r'\d+', test_string))))

    # print(dir_list)

    image_files = [os.path.join(image_folder,img) for img in dir_list if img.endswith(".png")]

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(os.path.join(movie_folder, movie_name), codec='mpeg4', audio=False, fps=fps)

def main():
    print(os.getcwd())
    image_folder = "fits_images"
    make_movies(image_folder=image_folder, movie_folder="./", movie_name="xrayproj.mp4")

if __name__ == "__main__":
    main()
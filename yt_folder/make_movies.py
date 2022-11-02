import os
import moviepy.video.io.ImageSequenceClip
import re
import numpy as np

def make_movies(image_folder, movie_name, **kwargs):
    image_folder
    fps=kwargs.get('fps', 10)
    suffix = kwargs.get('suffix', '.png')
    movie_folder = kwargs.get('movie_folder', image_folder)

    dir_list = os.listdir(image_folder)

    dir_list.sort(key=lambda test_string : list(
        map(int, re.findall(r'\d+', test_string))))


    image_files = [os.path.join(image_folder,img)
                for img in dir_list
                if img.endswith(suffix)]

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(os.path.join(movie_folder, movie_name), bitrate="12000k")
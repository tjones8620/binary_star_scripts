import os
import moviepy.video.io.ImageSequenceClip
image_folder='/mnt/local/thomas/wr140-cool-covertex/SimulationPlots_Thomas/WR140_hydro_cool_d3l6n064/Density/XY'
fps=10

image_files = [os.path.join(image_folder,img)
               for img in sorted(os.listdir(image_folder))
               if img.endswith(".png")]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('my_video.mp4')
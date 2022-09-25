from pathlib import Path
import os

home = str(Path.home())
print(home)

directory = "my/directory"

print(os.path.join(home, directory))

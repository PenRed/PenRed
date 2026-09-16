cd penred

# For Python 3.11 (Blender 4.2/4.5)
python3.11 -m bbext --all-wheels --split-platforms --no-cache --python 3.11
mv dist ../dist-311
rm -fr build

# For Python 3.13 (Blender 5.0+)
python3.13 -m bbext --all-wheels --split-platforms --no-cache --python 3.13
mv dist ../dist-313
rm -fr build

cd ..

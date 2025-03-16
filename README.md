# Conda-way of installing the required packages

1. Either use [direnv](https://direnv.net/) or manually activate the content of the `.envrc` file.

2. Create a conda/mamba environment (in all further commands, `mamba` can be replaced with `conda`, if you are using conda):

```bash
mamba create -n raisim python==3.10.12
mamba activate raisim
```

3. Install cxx compiler and related packages:

```bash
mamba install cxx-compiler ninja cmake eigen ffmpeg minizip
```

4. Install through cmake:

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DRAISIM_EXAMPLE=ON -DRAISIM_PY=ON -DPYTHON_EXECUTABLE=$(which python) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make install -j
```

By now raisim should be accessible and installation will be limited to your conda environment, thus, no need to worry about system-wide installation.

# Examples

You can run example from the `build` directory as executables. For example, to run `panda_example`:

```bash
cd build
./panda_example
```

# Visualizers

All the visualizers are available in the `third_party/raisimLib/{raisimUnity,raisimUnityOpengl}` directories. But you can start them with shortcats in `scripts` directory. The shortcuts also support differnt platforms. Note that `raisimUnityOpengl` is not supported by any Mac OS X, feel free to use `raisimUnity`.

```bash
bash scripts/raisimUnity.sh
bash scripts/raisimUnityOpengl.sh
```

## VSCode setup

VSCode can still be used for C++ project, but the suggestion is to use the following extensions instead of the default C++ ones:

## Known issues

1. While was tested on Mac OS X with Apple Sillicon, some examples failed to build with the same error:

```
error: no member named 'getSensorSet' in 'raisim::ArticulatedSystem'
```

Yet you're free to ignore it, snice it does not affect examples from the course.

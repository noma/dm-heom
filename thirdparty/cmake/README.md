This folder contains non-standard CMake additions.

- rpavlik contains a modified module for handling Git version information, originally from from https://github.com/rpavlik/cmake-modules
	- the used [fork](https://github.com/noma/cmake-modules) has an added `git_local_changes` function to detect local changes

```bash
cd rpavlik
wget https://github.com/noma/cmake-modules/raw/master/GetGitRevisionDescription.cmake
wget https://github.com/noma/cmake-modules/raw/master/GetGitRevisionDescription.cmake.in
wget https://github.com/noma/cmake-modules/raw/master/LICENSE_1_0.txt
cd ..
```

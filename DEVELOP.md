
# Release Notes

v1.1.0 - First PIP release, updated readme docs
v1.0.0 - Initial release of open source code base

# Instructions for developers

<hr>

## Testing

Several unit tests were written for py.test. From your virtual environment, run the following command from the repository root folder.

```bash
python -m pytest
```

Tests are automatically run via CircleCI. Results are located at [https://circleci.com/gh/synthego-open/ice](https://circleci.com/gh/synthego-open/ice).

<hr>

## Building the docker image locally

### Requirements:

* Docker http://docker.com

From a command line

```bash
git clone git:// # clone the repo
cd ice
docker build . -t ice # build the image from the Dockerfile
```

#### Running the examples

Running a single ICE analysis:

```bash
docker run -it -v ${PWD}:/data -w /ice -i ice:latest \
	python ice_analysis_single.py\
	--control /data/ice/tests/test_data/good_example_control.ab1 \
	--edited /data/ice/tests/test_data/good_example_edited.ab1 \
	--target AACCAGTTGCAGGCGCCCCA \
	--out /data/results/testing \
	--verbose
```

Running a batch analysis:

```bash
docker run -it -v ${PWD}:/data -w /ice -i ice:latest \
	python ice_analysis_batch.py\
	--in /data/ice/tests/test_data/batch_example.xlsx \
	--data /data/ice/tests/test_data/ \
	--out /data/results/ \
	--verbose
```

<hr>

## Installing locally from repository

#### Requirements:

* Python 3.X (last tested on 3.6, via anaconda)

#### Installation

Install your favorite python3 virtual environment (virtualenv, conda). We'll use conda for this example.

```bash
conda create --name ice_env python=3 # create a python3 virtual environment
source activate ice_env

git clone git:// # clone the repo
cd ice

pip install -r requirements.txt # install the python dependencies
```

#### Running the examples

Running a single ICE analysis:

```bash
./ice_analysis_single.py \
	--control ./ice/tests/test_data/good_example_control.ab1  \
	--edited ./ice/tests/test_data/good_example_edited.ab1 \
	--target AACCAGTTGCAGGCGCCCCA \
	--out results/testing \
	--verbose
```

Running a batch analysis:

```bash
./ice_analysis_batch.py \
	--in ./ice/tests/test_data/batch_example.xlsx \
	--out ./results/ \
	--data ./ice/tests/test_data/
	--verbose
```

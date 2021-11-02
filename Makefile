init:
	pip install --user -r requirements.txt

clean:
	rm -rf build
	rm -rf dist
	rm -rf overhang-assembly.egg-info
	rm -rf generate.egg-info

install:
	python setup.py install --user


develop:
	python setup.py develop --user

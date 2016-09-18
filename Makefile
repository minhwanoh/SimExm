clean:
	rm -rf ./build;
	find ./src/simulation/models -iname "*.c" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;

run:
	cd examples; python run_sim.py; cd ..

test:
	cd examples; python test.py; cd ..

build:
	python setup.py build_ext --inplace

install:
	make clean; make build

notebook:
	cd notebooks; jupyter notebook --ip='localhost' --port=8888; cd ..

clean:
	rm -rf ./build
	find ./database/methods -iname "*.c" -exec rm {} \;
	find ./simulation/models -iname "*.c" -exec rm {} \;
	find ./simulation -iname "methods.c" -exec rm {} \;
	find ./visualization -iname "methods.c" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;

run:
	cd scripts; python run_sim_easy.py; cd ..

runa:
	cd scripts; python run_sim_advanced.py; cd ..

build:
	python setup.py build_ext --inplace

install:
	make clean; make build

notebook:
	cd notebooks; jupyter notebook --ip='localhost' --port=8888; cd ..

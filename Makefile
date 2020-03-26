

cpp:
	./build.sh

python2:
	export SAGE_ROOT="${HOME}/sage/"; \
	export DOT_SAGE="${HOME}/.sage/"; \
	${SAGE_ROOT}/sage setup.py build_ext --inplace

python3:
	export SAGE_ROOT="${HOME}/sage3/sage/"; \
	export DOT_SAGE="${HOME}/sage3/.sage/"; \
	export SAGE_PYTHON3=yes; \
	${SAGE_ROOT}/sage setup.py build_ext --inplace

clean:
	\rm -f sagemath/*.cpp sagemath/*.c sagemath/*.html sagemath/*.so
	\rm -f *.cpp *.so *.html
	\rm -rf build/
	\rm -rf */__pycache__

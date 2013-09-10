NVCC=nvcc
NVCC_FLAGS=-Wall
CXX=g++
CXX_FLAGS=-Wall -std=c++0x

LIBS=-llog4cxx -lblas -llapack -latlas -lfftw3 -lfftw3f

src_dir=../src
run_dir=../run

sources=main config
bases_sources=collocation element messenger plan solver
one_d_sources=advection_one_d collocation_one_d diffusion_one_d element_one_d fftw_one_d solver_one_d
utils_sources=interpolate io solver_lapack utils_lapack

objects=$(sources:%=%.o)
	
bases_objects=$(bases_sources:%=bases/%.o)

one_d_objects=$(one_d_sources:%=one_d/%.o)

utils_objects=$(utils_sources:%=utils/%.o)

cuda_one_d_sources=diffusion_one_d_cuda element_one_d_cuda fftw_one_d_cuda solver_one_d_cuda
cuda_utils_sources=solver_utils_cuda utils_cublas

cuda_one_d_objects=$(cuda_one_d_sources:%=one_d/cuda/%.o)
cuda_utils_objects=$(cuda_utils_sources:%=utils/cuda/%.o)

main: $(objects) $(bases_objects) $(one_d_objects) $(utils_objects)
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(LIBS)
	
main_cuda: $(objects) $(bases_objects) $(one_d_objects) $(utils_objects) $(cuda_one_d_objects) $(cuda_utils_objects)
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(LIBS)

%.o: $(src_dir)/%.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

bases/%.o: $(src_dir)/bases/%.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

one_d/%.o: $(src_dir)/one_d/%.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

one_d/cuda/%.o: $(src_dir)/one_d/cuda/%.cu
	$(NVCC) $(NVCC_FLAGS) -o $@ -c $<

utils/%.o: $(src_dir)/utils/%.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

utils/cuda/%.o: $(src_dir)/utils/cuda/%.cu
	$(NVCC) $(NVCC_FLAGS) -o $@ -c $<

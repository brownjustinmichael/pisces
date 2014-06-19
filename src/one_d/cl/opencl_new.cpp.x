#include <iostream>
#include <../extern/cl.hpp>
#include <string>
#include <fstream>
#include <streambuf>
 
int main(){
    //get all platforms (drivers)
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    if(all_platforms.size()==0){
        std::cout<<" No platforms found. Check OpenCL installation!\n";
        exit(1);
    }
    cl::Platform default_platform=all_platforms[0];
    std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
 
    //get default device of the default platform
    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if(all_devices.size()==0){
        std::cout<<" No devices found. Check OpenCL installation!\n";
        exit(1);
    }
    cl::Device default_device=all_devices[1];
    std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
 
 
    cl::Context context({default_device});
    cl::Program::Sources sources;
	
	std::ifstream kernel_file ("../src/fft_utils.cl");
	std::string str ((std::istreambuf_iterator <char> (kernel_file)), std::istreambuf_iterator <char> ());
	sources.push_back ({str.c_str (), str.length ()});
	
    // kernel calculates for each element C=A+B
    cl::Program program (context, sources);
    if(program.build({default_device})!=CL_SUCCESS){
        std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
        exit(1);
    }
 
 
    // create buffers on the device
    cl::Buffer buffer_A(context,CL_MEM_READ_WRITE,sizeof(double)*10);
    cl::Buffer buffer_B(context,CL_MEM_READ_WRITE,sizeof(double)*20);
    cl::Buffer buffer_C(context,CL_MEM_READ_WRITE,sizeof(double)*20);
 
    double A[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double B[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1};
 
    //create queue to which we will push commands for the device.
    cl::CommandQueue queue(context,default_device);
 
    //write arrays A and B to the device
    queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0,sizeof(double)*10,A);
    queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0,sizeof(double)*10,B);
 
 
    //run the kernel
    cl::KernelFunctor simple_add(cl::Kernel(program,"simple_add"),queue,cl::NullRange,cl::NDRange(10),cl::NullRange);
    cl::KernelFunctor symmetrize(cl::Kernel(program,"symmetrize"),queue,cl::NullRange,cl::NDRange(10),cl::NullRange);
    // simple_add(buffer_A,buffer_B,buffer_C);
	symmetrize (10, buffer_A, buffer_C);
 
    //alternative way to run the kernel
    /*cl::Kernel kernel_add=cl::Kernel(program,"simple_add");
    kernel_add.setArg(0,buffer_A);
    kernel_add.setArg(1,buffer_B);
    kernel_add.setArg(2,buffer_C);
    queue.enqueueNDRangeKernel(kernel_add,cl::NullRange,cl::NDRange(10),cl::NullRange);
    queue.finish();*/
 
    double C[20];
    //read result C from the device to array C
    queue.enqueueReadBuffer(buffer_C,CL_TRUE,0,sizeof(double)*20,C);
 
    std::cout<<" result: \n";
    for(int i=0;i<20;i++){
        std::cout<<C[i]<<" ";
    }
	std::cout << "\n";
 
    return 0;
}
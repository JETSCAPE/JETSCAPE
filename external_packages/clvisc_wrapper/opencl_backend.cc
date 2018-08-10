/*******************************************************************************
 * Copyright (c) 2018-2019 LongGang Pang, lgpang@qq.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and/or associated documentation files (the
 * "Materials"), to deal in the Materials without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Materials, and to
 * permit persons to whom the Materials are furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Materials.
 *
 * THE MATERIALS ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * MATERIALS OR THE USE OR OTHER DEALINGS IN THE MATERIALS.
 ******************************************************************************/

#include "include/opencl_backend.h"

namespace clvisc {

CompileOption::CompileOption(){
    Define("USE_SINGLE_PRECISION");
};

CompileOption::CompileOption(bool use_single_precision, bool optimize) {
    if (use_single_precision) {
        Define("USE_SINGLE_PRECISION");
    }
    if (! optimize) {
        Define("-cl-opt-disable");
    }
}

void CompileOption::Define(std::string definition) {
    opt << "-D " << definition <<" ";
}

void CompileOption::KernelIncludePath(std::string abs_path) {
    opt << "-I " << abs_path <<" ";
}

void CompileOption::SetIntConst(std::string key, int value) {
    opt << "-D " << key << "=" << value << " ";
}

// for float values, use "#define key 0.33f" if value == 0.33.
void CompileOption::SetFloatConst(std::string key, float value) {
    opt << "-D " << key << "=" << std::setprecision(12) << std::fixed <<  value << "f ";
}

// for double values, use "#define key 0.33" if value == 0.33.
void CompileOption::SetDoubleConst(std::string key, double value) {
    opt << "-D " << key << "=" << value << " ";
}

std::string CompileOption::str() {
    return opt.str();
}

OpenclBackend::OpenclBackend(std::string device_type, int device_id) {
    // select device type and device id (if there are multiple cpu/gpus)
    if (device_type == "cpu" || device_type == "CPU") {
        device_type_ = CL_DEVICE_TYPE_CPU;
    } else if (device_type == "gpu" || device_type == "GPU") {
        device_type_ = CL_DEVICE_TYPE_GPU;
    } else {
        device_type_ = CL_DEVICE_TYPE_ALL;
    };
    device_id_ = device_id;
    // create context for the designated device type
    context_ = CreateContext_(device_type_);
    // choose one device if there are many of the same kind
    devices_ = context_.getInfo<CL_CONTEXT_DEVICES>();
    auto num_of_devices = devices_.size();
    if (device_id_ < 0 || device_id_ > num_of_devices-1) {
        DeviceInfo();
        throw std::out_of_range("device_id out of range");
    } else {
        device_ = devices_[device_id_];
    }
    queue_ = cl::CommandQueue(context_, device_, CL_QUEUE_PROFILING_ENABLE);
}



/** get the kernel excution time in units of seconds */
float OpenclBackend::ExcutionTime(cl::Event & event)
{
    cl_ulong tstart, tend;
    event.getProfilingInfo(CL_PROFILING_COMMAND_START, & tstart);
    event.getProfilingInfo(CL_PROFILING_COMMAND_END, & tend);
    //std::cout<<"#run time="<<(tend - tstart )/1000<<"ms\n";
    return (tend - tstart) * 1.0E-9 ;
}

cl::Context OpenclBackend::CreateContext_(const cl_int & device_type)
{
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.size() == 0) {
        std::cerr<<"No platform found, install CUDA or AMD SDK first\n";
        exit(-1);
    } else {
        for (int i=0; i < platforms.size(); i++) {
            std::vector<cl::Device> supportDevices;
            platforms.at(i).getDevices(CL_DEVICE_TYPE_ALL, &supportDevices);
            for (int j=0; j < supportDevices.size(); j++) {
                if (supportDevices.at(j).getInfo<CL_DEVICE_TYPE>() == device_type) {
                    //std::cout<<"#Found device "<<device_type<<" on platform "<<i<<std::endl;
                    cl_context_properties properties[] =
                    { CL_CONTEXT_PLATFORM, 
                        (cl_context_properties) (platforms.at(i))(),
                        0 };
                    return cl::Context(device_type, properties);
                }// Found supported device and platform
            }// End for devices
        }// End for platform
        //// if no platform support device_type, exit
        std::cerr<<"no platform support device type"<<device_type<<std::endl;
        exit(-1);
    }
}

cl::Program OpenclBackend::BuildProgram(std::string fname,
                                        const std::string & compile_option)
{ //// build programs and print the compile error if there is
    std::ifstream kernelFile(fname.c_str());
    if(!kernelFile.is_open()) {
        throw std::runtime_error("Fail to open kernel file: "+fname);
    }
    std::string sprog(std::istreambuf_iterator<char> (kernelFile),
                      (std::istreambuf_iterator<char> ()));
    cl::Program::Sources prog(1, std::make_pair(sprog.c_str(), sprog.length()));
    auto program = cl::Program(context_, prog);
    //programs.push(program);
    try{
        program.build(devices_, compile_option.c_str());
        kernelFile.close();
    } catch(cl::Error & err) {
        std::cerr << err.what() << "(" << err.err() << ")\n" \
            << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_);
    }

    return program;
}

cl::Buffer OpenclBackend::CreateBuffer(size_t bytes_of_buffer) {
    return cl::Buffer(context_, CL_MEM_READ_WRITE, bytes_of_buffer);
}

template <typename ValueType>
cl::Buffer OpenclBackend::CreateBufferByCopyVector(std::vector<ValueType> & source_vector,
                          bool read_only) {
    //copy content from a source vector to global memory of device
    if (read_only) {
        return cl::Buffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                source_vector.size()*sizeof(ValueType), source_vector.data());
    } else {
        return cl::Buffer(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                source_vector.size()*sizeof(ValueType), source_vector.data());
    }
}

cl::Image2D OpenclBackend::CreateImage2DByCopyVector(std::vector<cl_float4> & source_vector,
             size_t width, size_t height, bool read_only) {
    //copy content from a source vector to global memory of device
    cl::ImageFormat img_fmt;
    img_fmt.image_channel_order = CL_RGBA;
    img_fmt.image_channel_data_type = CL_FLOAT;
    size_t row_pitch = 0;
    if (read_only) {
        return cl::Image2D(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                img_fmt, width, height, row_pitch, source_vector.data());
    } else {
        return cl::Image2D(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                img_fmt, width, height, row_pitch, source_vector.data());
    }
}

void OpenclBackend::enqueue_run(const cl::Kernel & kernel_,
                                const cl::NDRange & global_size,
                                const cl::NDRange & local_size) {
    cl::Event event;
    queue_.enqueueNDRangeKernel(
            kernel_,                                // kernel name
            cl::NullRange,                          // offset 
            global_size,      // global size
            local_size,         // local size (automatically set by system)
            NULL,                     // event waitting list
            &event);       // event for profiling
    event.wait();
}


// from std::vector to cl::Buffer
template <typename ValueType>
void OpenclBackend::enqueue_copy(const std::vector<ValueType> & source_vector,
                                 cl::Buffer & dst_buffer)
{
    cl::Event event;
    queue_.enqueueWriteBuffer(
            dst_buffer,               // dst buffer
            CL_TRUE,                  // blocking reading
            0,                        // offset
            source_vector.size()*sizeof(ValueType),  // size
            source_vector.data(),     // source vector
            NULL,
            &event);                  
    event.wait();
}

// from cl::Buffer to std::vector
template <typename ValueType>
void OpenclBackend::enqueue_copy(const cl::Buffer & source_buffer,
                                 std::vector<ValueType> & dst_vector)
{
    cl::Event event;
    queue_.enqueueReadBuffer(
            source_buffer,            // source buffer
            CL_TRUE,                  // blocking reading
            0,                        // offset
            dst_vector.size()*sizeof(ValueType),  // size
            dst_vector.data(),       // dst vector
            NULL,
            &event);              
    event.wait();
}

// from cl::Buffer to cl::Buffer
void OpenclBackend::enqueue_copy(const cl::Buffer & source_buffer,
                                 cl::Buffer & dst_buffer,
                                 size_t size_in_bytes)
{
    cl::Event event;
    queue_.enqueueCopyBuffer(
            source_buffer,            // source buffer
            dst_buffer,               // dst buffer
            0,                        // src offset
            0,                        // dst offset
            size_in_bytes,            // size
            NULL,                     // waiting event-list
            &event);                  // event
    event.wait();
}


void OpenclBackend::DeviceInfo() {
    int device_id = 0;
    for (auto device : devices_) {
        std::cout << "Device ID: " << device_id << std::endl;
        std::cout << "Device Name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        std::cout << "Max computing units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
        std::cout << std::endl;
        std::cout << "Max workgroup size: " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
        std::cout << std::endl;
        std::cout << "Max work items in one work group: ";
        for (auto sz : device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()) {
            std::cout << sz << " ";
        }
        std::cout << std::endl;
        std::cout << "Global memory size: " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/1024/1024/1024 << "GB";
        std::cout << std::endl;
        std::cout << "Local memory size: " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/1024 << "KB";

        std::cout << std::endl << std::endl;
        device_id ++;
    }
}

// template member functions need explicit declearation on mac
template cl::Buffer OpenclBackend::CreateBufferByCopyVector(std::vector<cl_int> & source_vector, bool read_only);

template cl::Buffer OpenclBackend::CreateBufferByCopyVector(std::vector<cl_real> & source_vector, bool read_only);

template cl::Buffer OpenclBackend::CreateBufferByCopyVector(std::vector<cl_real4> & source_vector, bool read_only);

// cl_real3 is the same datatype as cl_real4 in cl.hpp, so one can not re-declear
//template cl::Buffer OpenclBackend::CreateBufferByCopyVector(std::vector<cl_real3> & source_vector, bool read_only);

template cl::Buffer OpenclBackend::CreateBufferByCopyVector(std::vector<cl_real8> & source_vector, bool read_only);

template void OpenclBackend::enqueue_copy(const std::vector<cl_int> & source_vector,  cl::Buffer & dst_buffer);

template void OpenclBackend::enqueue_copy(const std::vector<cl_real> & source_vector,  cl::Buffer & dst_buffer);

template void OpenclBackend::enqueue_copy(const std::vector<cl_real4> & source_vector,  cl::Buffer & dst_buffer);

template void OpenclBackend::enqueue_copy(const std::vector<cl_real8> & source_vector,  cl::Buffer & dst_buffer);

template void OpenclBackend::enqueue_copy(const cl::Buffer & source_buffer, std::vector<cl_int> & dst_vector);

template void OpenclBackend::enqueue_copy(const cl::Buffer & source_buffer, std::vector<cl_real> & dst_vector);

template void OpenclBackend::enqueue_copy(const cl::Buffer & source_buffer, std::vector<cl_real4> & dst_vector);

template void OpenclBackend::enqueue_copy(const cl::Buffer & source_buffer, std::vector<cl_real8> & dst_vector);

} // end namespace clvisc

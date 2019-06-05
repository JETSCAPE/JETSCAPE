#include"include/error_msgs.h"

std::string ErrorMessage(int err) {
    switch (err) {

        case 0:    return "CL_SUCCESS The sweet spot. ";

        case -1:    return "CL_DEVICE_NOT_FOUND clGetDeviceIDs if no OpenCL devices that matched device_type were found. ";

        case -2:    return "CL_DEVICE_NOT_AVAILABLE clCreateContext if a device in devices is currently not available even though "
                        "the device was returned by clGetDeviceIDs. ";

        case -3:    return "CL_COMPILER_NOT _AVAILABLE clBuildProgram if program is created with clCreateProgramWithSource and a "
                        "compiler is not available i.e. CL_DEVICE_COMPILER_AVAILABLE specified in the table of OpenCL Device Queries for "
                        "clGetDeviceInfo is set to CL_FALSE. ";

        case -4:    return "CL_MEM_OBJECT _ALLOCATION_FAILURE if there is a failure to allocate memory for buffer object. ";

        case -5:    return "CL_OUT_OF_RESOURCES if there is a failure to allocate resources required by the OpenCL implementation on "
                        "the device. ";

        case -6:    return "CL_OUT_OF_HOST_MEMORY if there is a failure to allocate resources required by the OpenCL implementation "
                        "on the host. ";

        case -7:    return "CL_PROFILING_INFO_NOT _AVAILABLE clGetEventProfilingInfo if the CL_QUEUE_PROFILING_ENABLE flag is "
                        "not set for the command-queue, if the execution status of the command identified by event is "
                        "not CL_COMPLETE or if event is a user event object. ";

        case -8:    return "CL_MEM_COPY_OVERLAP clEnqueueCopyBuffer, clEnqueueCopyBufferRect, clEnqueueCopyImage if "
                        "src_buffer and dst_buffer are the same buffer or subbuffer object and the source and destination regions "
                        "overlap or if src_buffer and dst_buffer are different sub-buffers of the same associated buffer object and "
                        "they overlap. The regions overlap if src_offset ≤ to dst_offset ≤ to src_offset + size – "
                        "1, or if dst_offset ≤ to src_offset ≤ to dst_offset + size – 1. ";

        case -9:    return "CL_IMAGE_FORMAT _MISMATCH clEnqueueCopyImage if src_image and dst_image do not use the same image "
                        "format. ";

        case -10:    return "CL_IMAGE_FORMAT_NOT _SUPPORTED clCreateImage if the image_format is not supported. ";

        case -11:    return "CL_BUILD_PROGRAM _FAILURE clBuildProgram if there is a failure to build the program executable. This "
                        "error will be returned if clBuildProgram does not return until the build has completed. ";

        case -12:    return "CL_MAP_FAILURE clEnqueueMapBuffer, clEnqueueMapImage if there is a failure to map the requested region "
                        "into the host address space. This error cannot occur for image objects created with CL_MEM_USE_HOST_PTR or "
                        "CL_MEM_ALLOC_HOST_PTR. ";

        case -13:    return "CL_MISALIGNED_SUB _BUFFER_OFFSET if a sub-buffer object is specified as the value for an argument that "
                        "is a buffer object and the offset specified when the sub-buffer object is created is not "
                        "aligned to CL_DEVICE_MEM_BASE_ADDR_ALIGN value for device associated with queue. ";

        case -14:    return "CL_EXEC_STATUS_ERROR_ FOR_EVENTS_IN_WAIT_LIST if the execution status of any of the events in "
                        "event_list is a negative integer value. ";

        case -15:    return "CL_COMPILE_PROGRAM _FAILURE clCompileProgram if there is a failure to compile the program source. This "
                        "error will be returned if clCompileProgram does not return until the compile has completed. ";

        case -16:    return "CL_LINKER_NOT_AVAILABLE clLinkProgram if a linker is not available i.e. CL_DEVICE_LINKER_AVAILABLE "
                        "specified in the table of allowed values for param_name for clGetDeviceInfo is set to CL_FALSE. ";

        case -17:    return "CL_LINK_PROGRAM_FAILURE clLinkProgram if there is a failure to link the compiled binaries and/or "
                        "libraries. ";

        case -18:    return "CL_DEVICE_PARTITION _FAILED clCreateSubDevices if the partition name is supported by the "
                        "implementation but in_device could not be further partitioned. ";

        case -19:    return "CL_KERNEL_ARG_INFO _NOT_AVAILABLE clGetKernelArgInfo if the argument information is not available "
                        "for kernel. ";

        case -30:    return "CL_INVALID_VALUE clGetDeviceIDs, clCreateContext This depends on the function: two or more coupled "
                        "parameters had errors. ";

        case -31:    return "CL_INVALID_DEVICE_TYPE clGetDeviceIDs if an invalid device_type is given ";

        case -32:    return "CL_INVALID_PLATFORM clGetDeviceIDs if an invalid platform was given ";

        case -33:    return "CL_INVALID_DEVICE clCreateContext, clBuildProgram if devices contains an invalid device or are not "
                        "associated with the specified platform. ";

        case -34:    return "CL_INVALID_CONTEXT if context is not a valid context. ";

        case -35:    return "CL_INVALID_QUEUE_PROPERTIES clCreateCommandQueue if specified command-queue-properties are valid "
                        "but are not supported by the device. ";

        case -36:    return "CL_INVALID_COMMAND_QUEUE if command_queue is not a valid command-queue. ";

        case -37:    return "CL_INVALID_HOST_PTR clCreateImage, clCreateBuffer This flag is valid only if host_ptr is not NULL. If "
                        "specified, it indicates that the application wants the OpenCL implementation to allocate memory for the memory "
                        "object and copy the data from memory referenced by host_ptr.CL_MEM_COPY_HOST_PTR and CL_MEM_USE_HOST_PTR are "
                        "mutually exclusive.CL_MEM_COPY_HOST_PTR can be used with CL_MEM_ALLOC_HOST_PTR to initialize the contents of "
                        "the cl_mem object allocated using host-accessible (e.g. PCIe) memory. ";

        case -38:    return "CL_INVALID_MEM_OBJECT if memobj is not a valid OpenCL memory object. ";

        case -39:    return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR if the OpenGL/DirectX texture internal format does not map to a "
                        "supported OpenCL image format. ";

        case -40:    return "CL_INVALID_IMAGE_SIZE if an image object is specified as an argument value and the image dimensions "
                        "(image width, height, specified or compute row and/or slice pitch) are not supported by device associated "
                        "with queue. ";

        case -41:    return "CL_INVALID_SAMPLER clGetSamplerInfo, clReleaseSampler, clRetainSampler, clSetKernelArg if "
                        "sampler is not a valid sampler object. ";

        case -42:    return "CL_INVALID_BINARY clCreateProgramWithBinary, clBuildProgram ";

        case -43:    return "CL_INVALID_BUILD_OPTIONS clBuildProgram if the build options specified by options are invalid. ";

        case -44:    return "CL_INVALID_PROGRAM if program is a not a valid program object. ";

        case -45:    return "CL_INVALID_PROGRAM_EXECUTABLE if there is no successfully built program executable available for "
                        "device associated with command_queue. ";

        case -46:    return "CL_INVALID_KERNEL_NAME clCreateKernel if kernel_name is not found in program. ";

        case -47:    return "CL_INVALID_KERNEL_DEFINITION clCreateKernel if the function definition for __kernel function given "
                        "by kernel_name such as the number of arguments, the argument types are not the same for "
                        "all devices for which the program executable has been built. ";

        case -48:    return "CL_INVALID_KERNEL if kernel is not a valid kernel object. ";

        case -49:    return "CL_INVALID_ARG_INDEX clSetKernelArg, clGetKernelArgInfo if arg_index is not a valid argument index. ";

        case -50:    return "CL_INVALID_ARG_VALUE clSetKernelArg, clGetKernelArgInfo if arg_value specified is not a valid value. ";

        case -51:    return "CL_INVALID_ARG_SIZE clSetKernelArg if arg_size does not match the size of the data type for an "
                        "argument that is not a memory object or if the argument is a memory object and "
                        "arg_size != sizeof(cl_mem) or if arg_size is zero and the argument is declared with the __local "
                        "qualifier or if the argument is a sampler and arg_size != sizeof(cl_sampler). ";

        case -52:    return "CL_INVALID_KERNEL_ARGS if the kernel argument values have not been specified. ";

        case -53:    return "CL_INVALID_WORK_DIMENSION if work_dim is not a valid value (i.e. a value between 1 and 3). ";

        case -54:    return "CL_INVALID_WORK_GROUP_SIZE if local_work_size is specified and number of work-items specified by "
                        "global_work_size is not evenly divisable by size of work-group given by local_work_size or does not match "
                        "the work-group size specified for kernel using the __attribute__ ((reqd_work_group_size(X, Y, Z))) "
                        "qualifier in program source.if local_work_size is specified and the total number of work-items in the work-group "
                        "computed as local_work_size[0] *… local_work_size[work_dim – 1] is greater than the value specified by "
                        "CL_DEVICE_MAX_WORK_GROUP_SIZE in the table of OpenCL Device Queries for clGetDeviceInfo.if local_work_size is NULL and the "
                        "__attribute__ ((reqd_work_group_size(X, Y, Z))) qualifier is used to declare the work-group size for kernel in the "
                        "program source. ";

        case -55:    return "CL_INVALID_WORK_ITEM_SIZE if the number of work-items specified in any of local_work_size[0], … "
                        "local_work_size[work_dim – 1] is greater than the corresponding values specified by CL_DEVICE_MAX_WORK_ITEM_SIZES[0], …. "
                        "CL_DEVICE_MAX_WORK_ITEM_SIZES[work_dim – 1]. ";

        case -56:    return "CL_INVALID_GLOBAL_OFFSET if the value specified in global_work_size + the corresponding values in "
                        "global_work_offset for any dimensions is greater than the sizeof(size_t) for the device on which the kernel "
                        "execution will be enqueued. ";

        case -57:    return "CL_INVALID_EVENT_WAIT_LIST if event_wait_list is NULL and num_events_in_wait_list > 0, or "
                        "event_wait_list is not NULL and num_events_in_wait_list is 0, or if event objects in event_wait_list are not "
                        "valid events. ";

        case -58:    return "CL_INVALID_EVENT if event objects specified in event_list are not valid event objects. ";

        case -59:    return "CL_INVALID_OPERATION if interoperability is specified by setting CL_CONTEXT_ADAPTER_D3D9_KHR, "
                        "CL_CONTEXT_ADAPTER_D3D9EX_KHR or CL_CONTEXT_ADAPTER_DXVA_KHR to a non-NULL value, and interoperability with another graphics API is "
                        "also specified. (only if the cl_khr_dx9_media_sharing extension is supported). ";

        case -60:    return "CL_INVALID_GL_OBJECT if texture is not a GL texture object whose type matches texture_target, if the "
                        "specified miplevel of texture is not defined, or if the width or height of the specified "
                        "miplevel is zero. ";

        case -61:    return "CL_INVALID_BUFFER_SIZE clCreateBuffer, clCreateSubBuffer if size is 0.Implementations may return "
                        "CL_INVALID_BUFFER_SIZE if size is greater than the CL_DEVICE_MAX_MEM_ALLOC_SIZE value specified in the table of allowed values "
                        "for param_name for clGetDeviceInfo for all devices in context. ";

        case -62:    return "CL_INVALID_MIP_LEVEL OpenGL-functions if miplevel is greater than zero and the OpenGL implementation "
                        "does not support creating from non-zero mipmap levels. ";

        case -63:    return "CL_INVALID_GLOBAL_WORK_SIZE if global_work_size is NULL, or if any of the values specified in "
                        "global_work_size[0], …global_work_size [work_dim – 1] are 0 or exceed the range given by the sizeof(size_t) for "
                        "the device on which the kernel execution will be enqueued. ";

        case -64:    return "CL_INVALID_PROPERTY clCreateContext Vague error, depends on the function ";

        case -65:    return "CL_INVALID_IMAGE_DESCRIPTOR clCreateImage if values specified in image_desc are not valid or if "
                        "image_desc is NULL. ";

        case -66:    return "CL_INVALID_COMPILER_OPTIONS clCompileProgram if the compiler options specified by options are "
                        "invalid. ";

        case -67:    return "CL_INVALID_LINKER_OPTIONS clLinkProgram if the linker options specified by options are invalid. ";

        case -68:    return "CL_INVALID_DEVICE_PARTITION_COUNT clCreateSubDevices if the partition name specified in properties "
                        "is CL_DEVICE_PARTITION_BY_COUNTS and the number of sub-devices requested exceeds "
                        "CL_DEVICE_PARTITION_MAX_SUB_DEVICES or the total number of compute units requested exceeds CL_DEVICE_PARTITION_MAX_COMPUTE_UNITS for "
                        "in_device, or the number of compute units requested for one or more sub-devices is less than "
                        "zero or the number of sub-devices requested exceeds CL_DEVICE_PARTITION_MAX_COMPUTE_UNITS for in_device. ";

        case -69:    return "CL_INVALID_PIPE_SIZE clCreatePipe if pipe_packet_size is 0 or the pipe_packet_size exceeds "
                        "CL_DEVICE_PIPE_MAX_PACKET_SIZE value for all devices in context or if pipe_max_packets is 0. ";

        case -70:    return "CL_INVALID_DEVICE_QUEUE clSetKernelArg when an argument is of type queue_t when it’s not a valid device "
                        "queue object. ";

        case -1000:    return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR clGetGLContextInfoKHR, clCreateContext CL and GL not on the "
                        "same device (only when using a GPU). ";

        case -1001:    return "CL_PLATFORM_NOT_FOUND_KHR clGetPlatform No valid ICDs found ";

        case -1002:    return "CL_INVALID_D3D10_DEVICE_KHR clCreateContext, clCreateContextFromType if the Direct3D 10 device "
                        "specified for interoperability is not compatible with the devices against which the context is to be "
                        "created. ";

        case -1003:    return "CL_INVALID_D3D10_RESOURCE_KHR clCreateFromD3D10BufferKHR, clCreateFromD3D10Texture2DKHR, "
                        "clCreateFromD3D10Texture3DKHR If the resource is not a Direct3D 10 buffer or texture object ";

        case -1004:    return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR clEnqueueAcquireD3D10ObjectsKHR If a mem_object is "
                        "already acquired by OpenCL ";

        case -1005:    return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR clEnqueueReleaseD3D10ObjectsKHR If a mem_object is not "
                        "acquired by OpenCL ";

        case -1006:    return "CL_INVALID_D3D11_DEVICE_KHR clCreateContext, clCreateContextFromType if the Direct3D 11 device "
                        "specified for interoperability is not compatible with the devices against which the context is to be "
                        "created. ";

        case -1007:    return "CL_INVALID_D3D11_RESOURCE_KHR clCreateFromD3D11BufferKHR, clCreateFromD3D11Texture2DKHR, "
                        "clCreateFromD3D11Texture3DKHR If the resource is not a Direct3D 11 buffer or texture object ";

        case -1008:    return "CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR clEnqueueAcquireD3D11ObjectsKHR If a mem_object is "
                        "already acquired by OpenCL ";

        case -1009:    return "CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR clEnqueueReleaseD3D11ObjectsKHR If a ‘mem_object’ is not "
                        "acquired by OpenCL ";

        case -1010:    return "CL_INVALID_D3D9_DEVICE_NV CL_INVALID_DX9_DEVICE_INTEL clCreateContext, "
                        "clCreateContextFromType If the Direct3D 9 device specified for interoperability is not compatible with the devices against "
                        "which the context is to be created ";

        case -1011:    return "CL_INVALID_D3D9_RESOURCE_NV CL_INVALID_DX9_RESOURCE_INTEL clCreateFromD3D9VertexBufferNV, "
                        "clCreateFromD3D9IndexBufferNV, clCreateFromD3D9SurfaceNV, clCreateFromD3D9TextureNV, clCreateFromD3D9CubeTextureNV, "
                        "clCreateFromD3D9VolumeTextureNV If a ‘mem_object’ is not a Direct3D 9 resource of the required type ";

        case -1012:    return "CL_D3D9_RESOURCE_ALREADY_ACQUIRED_NV CL_DX9_RESOURCE_ALREADY_ACQUIRED_INTEL "
                        "clEnqueueAcquireD3D9ObjectsNV If any of the ‘mem_objects’ is currently already acquired by OpenCL ";

        case -1013:    return "CL_D3D9_RESOURCE_NOT_ACQUIRED_NV CL_DX9_RESOURCE_NOT_ACQUIRED_INTEL "
                        "clEnqueueReleaseD3D9ObjectsNV If any of the ‘mem_objects’ is currently not acquired by OpenCL ";

        case -1092:    return "CL_EGL_RESOURCE_NOT_ACQUIRED_KHR clEnqueueReleaseEGLObjectsKHR If a ‘mem_object’ is not acquired "
                        "by OpenCL ";

        case -1093:    return "CL_INVALID_EGL_OBJECT_KHR clCreateFromEGLImageKHR, clEnqueueAcquireEGLObjectsKHR If a "
                        "‘mem_object’ is not a EGL resource of the required type ";

        case -1094:    return "CL_INVALID_ACCELERATOR_INTEL clSetKernelArg when ‘arg_value’ is not a valid accelerator object, and "
                        "by clRetainAccelerator, clReleaseAccelerator, and clGetAcceleratorInfo when ‘accelerator’ is not a "
                        "valid accelerator object ";

        case -1095:    return "CL_INVALID_ACCELERATOR_TYPE_INTEL clSetKernelArg, clCreateAccelerator when ‘arg_value’ is not an "
                        "accelerator object of the correct type, or when ‘accelerator_type’ is not a valid accelerator type ";

        case -1096:    return "CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL clCreateAccelerator when values described by "
                        "‘descriptor’ are not valid, or if a combination of values is not valid ";

        case -1097:    return "CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL clCreateAccelerator when ‘accelerator_type’ is a valid "
                        "accelerator type, but it not supported by any device in ‘context’ ";

        case -1098:    return "CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL clCreateContext, clCreateContextFromType If the VA API "
                        "display specified for interoperability is not compatible with the devices against which the context is to "
                        "be created ";

        case -1099:    return "CL_INVALID_VA_API_MEDIA_SURFACE_INTEL clEnqueueReleaseVA_APIMediaSurfacesINTEL If ‘surface’ is "
                        "not a VA API surface of the required type, by clGetMemObjectInfo when ‘param_name’ is "
                        "CL_MEM_VA_API_MEDIA_SURFACE_INTEL when was not created from a VA API surface, and from clGetImageInfo when ‘param_name’ is "
                        "CL_IMAGE_VA_API_PLANE_INTEL and ‘image’ was not created from a VA API surface ";

        case -1100:    return "CL_VA_API_MEDIA_SURFACE_ALREADY_ACQUIRED_INTEL clEnqueueReleaseVA_APIMediaSurfacesINTEL If "
                        "any of the ‘mem_objects’ is already acquired by OpenCL ";

        case -1101:    return "CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL clEnqueueReleaseVA_APIMediaSurfacesINTEL If any of "
                        "the ‘mem_objects’ are not currently acquired by OpenCL ";

        case -9999:    return "NVidia clEnqueueNDRangeKernel Illegal read or write to a buffer ";

        default: return "Unknown err_id";
    }
}
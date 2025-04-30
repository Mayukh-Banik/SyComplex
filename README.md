syComplex

A temporary sycl capable Complex data type. ICPX currently supports std::complex, but I've also used HIPSYCL which doesn't so this alleviates this.

Functionality aimed to be completed:
https://github.khronos.org/SYCL_Reference/built-in-functions.html

```
https://github.khronos.org/SYCL_Reference/iface/function-objects.html
sycl::plus - done
sycl::multiplies - done
sycl::bit_and - done
sycl::bit_or - done
sycl::bit_xor - done
sycl::logical_and - done
sycl::logical_or - done
sycl::minimum - done
sycl::maximum - done

https://github.khronos.org/SYCL_Reference/iface/group-functions.html
No functionality from here will be implemented

https://github.khronos.org/SYCL_Reference/iface/group-algorithms-library.html
No functionatilty from here will be implemented

https://github.khronos.org/SYCL_Reference/iface/math-functions.html
sycl::acos - done
sycl::acosh - done
sycl::acospi - done
sycl::asin - done
sycl::asinh - done
sycl::asinpi - done
sycl::atan - done
sycl::atan2 - done
sycl::atanh - done
sycl::atanpi - done
sycl::atan2pi - done
sycl::cbrt - done
sycl::ceil - will not be implemented
sycl::copysign - done 


```
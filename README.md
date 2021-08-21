# ASM-Math
<p>Mathematical floating-point functions coded in x64 Assembly for NASM using the SSE instruction set.</p>
<p>ASM-Math aims to find a reasonable balance between accuracy and speed. This means
that all results should be correct. However, rounding errors might cause the last digit
to be wrong. The library also avoids using the stack/cache, instead it saves all variables in registers.</p>
<p>Currently there is no version of this library that uses the AVX instruction set and the x87 FPU.</p>

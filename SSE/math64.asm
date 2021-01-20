; math64.asm provides several scalar  mathematical functions in additions to the ones that come with
; the SSE and the "regular" x86 assembly instructions.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Copyright Johannes Kloimböck 2021.                                                            ;;
;; Distributed under the Boost Software License, Version 1.0.                                    ;;
;; (See accompanying file LICENSE or copy at                                                     ;;
;; https://www.boost.org/LICENSE_1_0.txt)                                                        ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

section .data

section .bss

section .text
%macro expSmall64 1	    ; calculates exp(x) where x is between -1.0 and 1.0 for 64-Bit float values
			            ; modifies the following registers: rdi, xmm1, xmm2, xmm3, xmm4, xmm5
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi	    ; move 1.0 into xmm1 (sum)
	movsd xmm2, xmm1	; term
	movsd xmm3, xmm2	; divisor for term
	movsd xmm4, xmm3	; used to add 1.0 to xmm3
	mov rdi, 0xbca0000000000000
	movq xmm5, rdi		; move negative term limit into xmm5
	mov rdi, 0x3ca0000000000000
	movq xmm6, rdi		; move positive term limit into xmm6
	
	.expLoop:		    ; calculating the Taylor series
		mulsd xmm2, %1
		divsd xmm2, xmm3
		addsd xmm3, xmm4
		addsd xmm1, xmm2
		ucomisd xmm2, xmm6
		jnc .expLoop
		ucomisd xmm2, xmm5
		jc .expLoop
%endmacro

exp64:			        ; calculates exp(x) for any double-precision floating-point input
			            ; modifies the following registers: rdi, rsi, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	
	cvtsd2si rdi, xmm0	; integer part of the exponent is now in rdi
	cvtsi2sd xmm4, rdi	; integer part of the exponent is also in xmm4
	movsd xmm2, xmm0	; store xmm0 in xmm2
	subsd xmm2, xmm4	; extract the fractional part of the exponent in xmm2
	
	mov rsi, 0x4005bf0a8b145769	; bits of e as a float in rsi
	movq xmm0, rsi			; e in xmm0
	
	call powfi64			; exp(rdi)
	
	movsd xmm7, xmm0		; store exp(rdi) in xmm7
	movsd xmm0, xmm2		; store the fractional exponent in xmm0
	
	expSmall64 xmm0			; exp(fractional Part) by using a macro
	movsd xmm0, xmm1
	mulsd xmm0, xmm7		; exp(int) * exp(fraction) in xmm0
	
	mov rsp, rbp
	pop rbp
	ret

powfi64:		            ; calculates x^n, where x is a real (double-precision) number (xmm0) and n is a 64-bit signed integer (rdi)
			                ; modifies the following registers: rdi, rsi, rcx, rdx, r8, r9, xmm0, xmm1
	mov rdx, 0x3ff0000000000000
	movsd xmm1, xmm0
	movq xmm0, rdx
	
	mov r8, 0x8000000000000000
	mov r9, rdi
	neg r9
	mov rsi, rdi
	test rdi, r8
	cmovnz rdi, r9
	mov r8, 0b1
	bsr rcx, rdi
	jz .zeroExponent
	inc rcx
	mov r9, 0b111111
	
	.powLoop:
		dec rcx
		test rdi, r8
		jz .zeroBit
		mulsd xmm0, xmm1
		.zeroBit:
			mulsd xmm1, xmm1
		shl r8, 0b1
		test rcx, r9
		jnz .powLoop
		
	mov r8, 0x8000000000000000
	test rsi, r8
	jz .zeroExponent
	movq xmm1, rdx
	divsd xmm1, xmm0
	movsd xmm0, xmm0
	
	.zeroExponent:

	ret 

%macro ldInprecise64 1		; estimates the size of ld(xmm0); only modifies rax for the result
	movq rax, %1
	shr rax, 52
	sub rax, 0x3ff
%endmacro

log64:				        ; calculates the natural logarithm of a double-precision floating-point number
				            ; returns NaN, if the argument is negative or NaN
				            ; returns negative infinity, if the argument is +0.0
				            ; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8 - r11, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	mov rdi, 0x8000000000000000
	movq rsi, xmm0
	test rsi, rdi
	jnz .nan
	xor rdi, rdi
	xor rsi, rdi
	jz .neginf
	
	ldInprecise64 xmm0	    ; use macro
	
	cvtsi2sd xmm1, rax
	mov rax, 0x3fe62e42fefa39ef	; log(2) in xmm2
	movq xmm2, rax
	mulsd xmm1, xmm2
	
	movq rax, xmm0

	mov r10, 3
	
	.newtonLoop:
		movsd xmm0, xmm1
		movq r11, xmm1
		call exp64
		movq xmm2, rax
		movsd xmm3, xmm2
		subsd xmm2, xmm0
		addsd xmm3, xmm0
		divsd xmm2, xmm3
		addsd xmm2, xmm2
		movq xmm1, r11
		addsd xmm1, xmm2
		dec r10
		jnz .newtonLoop

	movsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	
	ret
	
	.nan:
		mov rax, 0x7fffffffffffffff
		movq xmm0, rax
		mov rsp, rbp
		pop rbp
		ret
	
	.neginf:
		mov rax, 0xfff0000000000000
		movq xmm0, rax
		mov rsp, rbp
		pop rbp
		ret

log1064:		            ; calculates the common logarithm (base 10) of a double-precision floating-point number in xmm0
			                ; returns NaN, if the argument is negative or NaN
			                ; returns negative infinity, if the argument is +0.0
			                ; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8, r9, r10, r11, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	call log64
	mov rdi, 0x3fdbcb7b1526e50e	; log10(e)
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

log264:			            ; calculates the binary logarithm (base 2) of a double-precision floating-point number in xmm0
			                ; returns NaN, if the argument is negative or NaN
			                ; returns negative infinity, if the argument is +0.0
			                ; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8, r9, r10, r11, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	call log64
	mov rdi, 0x3ff71547652b82fe	; log2(e)
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

powff64:		            ; calculates x^y, where x (xmm0) and y (xmm1) are both double-precision floating-point numbers
			                ; returns NaN, if one argument is NaN, or xmm0 is negative and xmm1 is not an integer
			                ; modifies the following registers: rdi, (rsi, rax, rcx, rdx, r8, r9, r10, r11), xmm0, xmm1, (xmm2 - xmm7)
	push rbp
	mov rbp, rsp
	cvtsd2si rdi, xmm1
	cvtsi2sd xmm2, rdi
	ucomisd xmm2, xmm1
	jz .intPower
	movq rdi, xmm1
	push rdi
	call log64
	movsd xmm1, xmm0
	movq xmm0, rdi
	pop rdi
	mulsd xmm0, xmm1
	call exp64
	mov rsp, rbp
	pop rbp
	ret
	
	.intPower:
		call powfi64
		mov rsp, rbp
		pop rbp
		ret

sinh64:			            ; calculates the hyperbolic sine of a double-precision floating-point number
			                ; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8, r9, r10, r11, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	call exp64
	movsd xmm1, xmm0
	mov rdi, 0x3ff0000000000000
	movq xmm2, rdi
	divsd xmm2, xmm1
	subsd xmm0, xmm2
	mov rdi, 0x3fe0000000000000
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

cosh64:			            ; calculates the hyperbolic cosine of a double-precision floating-point number
			                ; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8, r9, r10, r11, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	call exp64
	movsd xmm1, xmm0
	mov rdi, 0x3ff0000000000000
	movq xmm2, rdi
	divsd xmm2, xmm1
	addsd xmm0, xmm2
	mov rdi, 0x3fe0000000000000
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

tanh64:			            ; calculates the hyperbolic tangent of a double-precision floating-point number
			                ; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8, r9, r10, r11, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	addsd xmm0, xmm0
	call exp64
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi
	movsd xmm2, xmm1
	addsd xmm1, xmm0
	subsd xmm0, xmm2
	divsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret
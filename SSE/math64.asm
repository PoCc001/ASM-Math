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
cbrt64:				; calculates the cube root of a double-precision floating-point number
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, xmm0 - xmm4
	movq rax, xmm0
	mov rsi, 0x8000000000000000
	and rsi, rax		; sign of the input is in rsi
	xor rax, rsi		; absolute value of the input in rax
	movq xmm0, rax
	xor rdi, rdi
	xor rdx, rdx
	xor rax, rdi		; check if the input is 0.0 or -0.0
	jz .zero
	mov rdi, 0x3
	div rdi			; divide by 3
	mov rdi, 0x2aa0000000000000
	mov rcx, 0x7ff0000000000000
	mov dl, 0x20
	and rcx, rax
	mov cl, 0x6
	cmovz cx, dx		; do more iterations, if the input is denormal
	add rax, rdi
	movq xmm1, rax
	mov rdi, 0x10000000000000
	mov rdx, 0x3fd5555555555555
	movq xmm3, rdx
	
	.newton:
		movq rdx, xmm1
		add rdx, rdi
		movq xmm2, rdx
		mulsd xmm1, xmm1
		movsd xmm4, xmm0
		divsd xmm4, xmm1
		addsd xmm2, xmm4
		mulsd xmm2, xmm3
		movsd xmm1, xmm2
		dec cl
		jnz .newton
	
	movq rdi, xmm1
	xor rdi, rsi
	movq xmm0, rdi
	
	ret
	
	.zero:
		movq rdi, xmm0
		xor rdi, rsi
		movq xmm0, rdi	
	ret

%macro expSmall64 1	    	; calculates exp(x) where x is between -1.0 and 1.0 for 64-Bit float values
				; modifies the following registers: rdi, xmm1, xmm2, xmm3, xmm4, xmm5
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi	    	; move 1.0 into xmm1 (sum)
	movsd xmm2, xmm1	; term
	movsd xmm3, xmm2	; divisor for term
	movsd xmm4, xmm3	; used to add 1.0 to xmm3
	xorpd xmm5, xmm5
	.expLoop:		; calculating the Taylor series
		mulsd xmm2, %1
		divsd xmm2, xmm3
		addsd xmm3, xmm4
		addsd xmm1, xmm2
		ucomisd xmm1, xmm5	; stop, when the sum doesn't change anymore
		movsd xmm5, xmm1
		jnz .expLoop
%endmacro

exp64:			        ; calculates exp(x) for any double-precision floating-point input
				; modifies the following registers: rdi, rsi, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	
	cvtsd2si rdi, xmm0	; integer part of the exponent is now in rdi
	cvtsi2sd xmm4, rdi	; integer part of the exponent is also in xmm4
	movsd xmm2, xmm0	; store xmm0 in xmm2
	subsd xmm2, xmm4	; extract the fractional part of the exponent in xmm2
	
	mov rsi, 0x4005bf0a8b145769	; bits of e as a float in rsi
	movq xmm0, rsi			; e in xmm0
	
	call powfi64			; exp(rdi)
	
	movsd xmm6, xmm0		; store exp(rdi) in xmm7
	movsd xmm0, xmm2		; store the fractional exponent in xmm0
	
	expSmall64 xmm0			; exp(fractional Part) by using a macro
	movsd xmm0, xmm1
	mulsd xmm0, xmm6		; exp(int) * exp(fraction) in xmm0
	
	mov rsp, rbp
	pop rbp
	ret
	
exp1m64:			; calculates exp(x) - 1 with x as a double-precision floating-point number
				; modifies the following registers: rdi, rsi, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	
	mov rdi, 0x7fffffffffffffff
	movq rsi, xmm0
	and rsi, rdi
	movq xmm1, rsi
	mov rdi, 0x3fd0000000000000
	movq xmm2, rdi
	ucomisd xmm1, xmm2
	jnc .useExp64
	
	movsd xmm1, xmm0
	movsd xmm3, xmm1
	mov rdi, 0x3ff0000000000000
	movq xmm2, rdi
	movsd xmm4, xmm2
	xorpd xmm5, xmm5
	
	.taylorLoop:
		mulsd xmm3, xmm1
		addsd xmm4, xmm2
		divsd xmm3, xmm4
		addsd xmm0, xmm3
		ucomisd xmm0, xmm5
		movsd xmm5, xmm0
		jnz .taylorLoop
	
	mov rsp, rbp
	pop rbp
	
	ret
	
	.useExp64:
		call exp64
		mov rdi, 0xbff0000000000000
		movq xmm1, rdi
		addsd xmm0, xmm1
		
	mov rsp, rbp
	pop rbp
	
	ret

powfi64:			; calculates x^n, where x is a real (double-precision) number (xmm0) and n is a 64-bit signed integer (rdi)
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
	mov r8d, 0b1
	bsr rcx, rdi
	jz .zeroExponent
	inc ecx
	mov r9d, 0b111111
	
	.powLoop:
		sub ecx, 0b1
		test rdi, r8
		jz .zeroBit
		mulsd xmm0, xmm1
		.zeroBit:
			mulsd xmm1, xmm1
		shl r8, 0b1
		test ecx, r9d
		jnz .powLoop
		
	mov r8, 0x8000000000000000
	test rsi, r8
	jz .zeroExponent
	movq xmm1, rdx
	divsd xmm1, xmm0
	movsd xmm0, xmm0
	
	.zeroExponent:

	ret 

log64:				; calculates the natural logarithm of a double-precision floating-point number
				; returns NaN, if the argument is negative or NaN
				; returns negative infinity, if the argument is +0.0
				; modifies the following registers: rdi, rsi, rcx, rdx, xmm0 - xmm6
	mov rdi, 0x8000000000000000
	movq rsi, xmm0
	and rdi, rsi
	jnz .nan
	xor edi, edi
	xor rdi, rsi
	jz .neginf
	
	movq rdx, xmm0
	shr rdx, 52
	jnz .notDenormal
	
	movq rsi, xmm0
	bsr rdi, rsi
	sub edi, 0b110101
	not edi
	mov ecx, edi
	sub edx, edi
	
	shl rsi, cl
	movq xmm0, rsi
	inc edx
	
	.notDenormal:
		sub edx, 0x3ff
		movq rdi, xmm0
		mov rsi, 0x3fffffffffffffff
		and rdi, rsi
		mov rsi, 0x3ff0000000000000
		or rdi, rsi
		movq xmm0, rdi		; argument is now between 1 and 2
		sqrtsd xmm0, xmm0	; argument is now between 1 and sqrt(2)
	
	movq xmm1, rsi
	movsd xmm2, xmm1
	addsd xmm1, xmm0	; (x + 1)
	subsd xmm0, xmm2	; (x - 1)
	
	divsd xmm0, xmm1
	movsd xmm1, xmm0
	movsd xmm4, xmm0
	mulsd xmm4, xmm4	; (x - 1)^2 / (x + 1)^2
	
	movq xmm3, rsi		; factor
	movsd xmm2, xmm3
	addsd xmm2, xmm2	; incrementor
	
	mov rsi, 0xbff0000000000000
	movq xmm6, rsi
	
	.taylorLoop:
		mulsd xmm1, xmm4
		addsd xmm3, xmm2
		movsd xmm5, xmm1
		divsd xmm5, xmm3
		addsd xmm0, xmm5
		ucomisd xmm0, xmm6
		movsd xmm6, xmm0
		jnz .taylorLoop
	
	mov rdi, 0x20000000000000
	movq rsi, xmm0
	add rsi, rdi
	movq xmm0, rsi
	
	mov rdi, 0x3fe62e42fefa39ef
	movq xmm1, rdi
	movsx rdx, edx
	cvtsi2sd xmm2, rdx
	mulsd xmm1, xmm2
	addsd xmm0, xmm1
	
	ret
	
	.nan:
		mov rdi, 0x7fffffffffffffff
		movq xmm0, rdi
		
		ret
	
	.neginf:
		mov rdi, 0xfff0000000000000
		movq xmm0, rdi
		
		ret


log1064:			; calculates the common logarithm (base 10) of a double-precision floating-point number in xmm0
				; returns NaN, if the argument is negative or NaN
				; returns negative infinity, if the argument is +0.0
				; modifies the following registers: rdi, rsi, rcx, rdx, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	call log64
	mov rdi, 0x3fdbcb7b1526e50e	; log10(e)
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

log264:				; calculates the binary logarithm (base 2) of a double-precision floating-point number in xmm0
				; returns NaN, if the argument is negative or NaN
				; returns negative infinity, if the argument is +0.0
				; modifies the following registers: rdi, rsi, rcx, rdx, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	call log64
	mov rdi, 0x3ff71547652b82fe	; log2(e)
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

logb64:				; calculates the logarithm of x (xmm0) with respect to base y (xmm1) - x and y being two double-precision floating-point numbers
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	
	movsd xmm1, xmm7
	call log64
	movq rax, xmm0
	movsd xmm0, xmm7
	call log64
	movq xmm1, rax
	divsd xmm0, xmm1
	
	mov rsp, rbp
	pop rbp
	
	ret
	
log1p64:			; calculates log(1 + x), with x as a double-precision floating-point number.
				; modifies the following registers: rdi, rsi, rcx, rdx, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	
	mov rdi, 0x7fffffffffffffff
	movq rsi, xmm0
	and rsi, rdi
	movq xmm1, rsi
	mov rdi, 0x3fd0000000000000
	movq xmm2, rdi
	ucomisd xmm1, xmm2
	jnc .useLog64
	mov rdi, 0xbfefffffffffffff
	movq xmm1, rdi
	ucomisd xmm0, xmm1
	jc .useLog64
	
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi
	movsd xmm2, xmm1
	movsd xmm3, xmm0
	movsd xmm4, xmm3
	movq rdi, xmm3
	mov rsi, 0x8000000000000000
	or rdi, rsi
	movq xmm3, rdi
	mov rdi, 0x4000000000000000
	movq xmm6, rdi
	
	.taylorLoop:
		mulsd xmm4, xmm3
		movsd xmm5, xmm4
		addsd xmm1, xmm2
		divsd xmm5, xmm1
		addsd xmm0, xmm5
		ucomisd xmm0, xmm6
		movsd xmm6, xmm0
		jnz .taylorLoop
	
	mov rsp, rbp
	pop rbp
	
	ret
	
	.useLog64:
		mov rdi, 0x3ff0000000000000
		movq xmm1, rdi
		addsd xmm0, xmm1
		call log64
		
	mov rsp, rbp
	pop rbp
	
	ret

powff64:			; calculates x^y, where x (xmm0) and y (xmm1) are both double-precision floating-point numbers
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

sinh64:				; calculates the hyperbolic sine of a double-precision floating-point number
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

cosh64:				; calculates the hyperbolic cosine of a double-precision floating-point number
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

tanh64:				; calculates the hyperbolic tangent of a double-precision floating-point number
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

%macro sinTaylor64 1		; calculates the sine of a double-precision floating-point number using the Taylor series for the sine function
				; modifies the following registers: rdi, (xmm0), xmm1 - xmm6
	movsd xmm1, %1		; sum
	movsd xmm4, xmm1	; term
	mulsd %1, %1		; x*x
	xorpd xmm2, xmm2
	subsd xmm2, %1
	movsd %1, xmm2		; -x*x
	mov rdi, 0x3ff0000000000000
	movq xmm2, rdi
	movsd xmm3, xmm2	; 1.0 in xmm2, xmm3 and xmm5
	movsd xmm5, xmm3	; xmm3 is the incrementor
	mov rdi, 0xc000000000000000
	movq xmm6, rdi		; -2.0 in xmm6
	
	.taylorLoop:
		addsd xmm2, xmm3	; increment xmm2
		movsd xmm5, xmm2	; set up denominator of term
		addsd xmm2, xmm3	; increment xmm2 again
		mulsd xmm5, xmm2	; multiply xmm2 and xmm5 together
		divsd xmm4, xmm5	; divide xmm5 from the term
		mulsd xmm4, %1
		addsd xmm1, xmm4
		ucomisd xmm1, xmm6
		movsd xmm6, xmm1
		jnz .taylorLoop
	movsd xmm0, xmm1
%endmacro

sin64:				; calculates the sine of a double-precision floating-point number
				; modifies the following registers: rdi, rsi, xmm0 - xmm6
	mov rdi, 0x401921fb54442d18
	movq xmm1, rdi	; 2*pi in xmm1
	movsd xmm2, xmm0	; copy input value
	divsd xmm2, xmm1	; see how often 2 * pi in the input
	cvtsd2si rdi, xmm2	; make it an integer in rdi
	mov rcx, 0x43b0000000000000	; threshold
	movq xmm3, rcx		; threshold in xmm3
	movq rsi, xmm2		; store bits of xmm2 in rsi
	ucomisd xmm2, xmm3	; is xmm2 smaller than 2^60?
	cmovc rsi, rdi		; if so, copy the integer quotient into rsi
	cvtsi2sd xmm2, rsi	; quotient to xmm2
	mulsd xmm1, xmm2	; multiply the quotient with 2 * pi
	subsd xmm0, xmm1	; subtract it from the input
	
	sinTaylor64 xmm0
	ret

%macro cosTaylor64 1		; calculates the cosine of a double-precision floating-point number using the Taylor series for the cosine function
				; modifies the following registers: rdi, (xmm0), xmm1 - xmm6
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi		; sum of the series
	mulsd %1, %1		; x*x
	xorpd xmm2, xmm2
	subsd xmm2, %1
	movsd %1, xmm2		; -x*x
	movsd xmm3, xmm1	; factor for factorial in denominator
	movsd xmm4, xmm3	; incrementor of factor
	movsd xmm2, xmm4	; term
	mov rdi, 0xc000000000000000
	movq xmm5, rdi		; sum to test against
	movsd xmm6, xmm2	; single factor for factorial

	.taylorLoop:
		movsd xmm3, xmm6
		addsd xmm6, xmm4
		mulsd xmm3, xmm6
		addsd xmm6, xmm4
		divsd xmm2, xmm3
		mulsd xmm2, %1
		addsd xmm1, xmm2
		ucomisd xmm1, xmm5
		movsd xmm5, xmm1
		jnz .taylorLoop
	movsd xmm0, xmm1
%endmacro

cos64:				; calculates the cosine of a double-precision floating-point number
				; modifies the following registers: rdi, rsi, xmm0 - xmm6
	mov rdi, 0x401921fb54442d18
	movq xmm1, rdi	; 2*pi in xmm1
	movsd xmm2, xmm0	; copy input value
	divsd xmm2, xmm1	; see how often 2 * pi in the input
	cvtsd2si rdi, xmm2	; make it an integer in rdi
	mov rcx, 0x43b0000000000000	; threshold
	movq xmm3, rcx		; threshold in xmm3
	movq rsi, xmm2		; store bits of xmm2 in rsi
	ucomisd xmm2, xmm3	; is xmm2 smaller than 2^60?
	cmovc rsi, rdi		; if so, copy the integer quotient into rsi
	cvtsi2sd xmm2, rsi	; quotient to xmm2
	mulsd xmm1, xmm2	; multiply the quotient with 2 * pi
	subsd xmm0, xmm1	; subtract it from the input
	
	cosTaylor64 xmm0
	ret

tan64:				; calculates the tangent of a double-precision floating-point number
				; modifies the following registers: rdi, r8, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	
	movsd xmm7, xmm0
	call sin64
	movq r8, xmm0
	movsd xmm0, xmm7
	call cos64
	movq xmm1, r8
	divsd xmm1, xmm0	; tan(x) = sin(x) / cos(x)
	movsd xmm0, xmm1
	
	mov rsp, rbp
	pop rbp
	ret

arsinh64:			; calculates the inverse hyperbolic sine of a 64-bit floating-point number
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8 - r10, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	movsd xmm1, xmm0
	mulsd xmm1, xmm1
	mov rdi, 0x3ff0000000000000
	movq xmm2, rdi
	addsd xmm1, xmm2
	sqrtsd xmm1, xmm1
	addsd xmm0, xmm1
	call log64
	mov rsp, rbp
	pop rbp
	ret

arcosh64:			; calculates the inverse hyperbolic cosine of a 64-bit floating-point number
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8 - r10, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	movsd xmm1, xmm0
	mulsd xmm1, xmm1
	mov rdi, 0x3ff0000000000000
	movq xmm2, rdi
	subsd xmm1, xmm2
	sqrtsd xmm1, xmm1
	addsd xmm0, xmm1
	call log64
	mov rsp, rbp
	pop rbp
	ret

artanh64:			; calculates the inverse hyperbolic tangent of a 64-bit floating-point number
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, r8 - r10, xmm0 - xmm7
	push rbp
	mov rbp, rsp
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi
	addsd xmm1, xmm0
	movq xmm2, rdi
	subsd xmm2, xmm0
	divsd xmm1, xmm2
	movsd xmm0, xmm1
	call log64
	mov rdi, 0x3fe0000000000000
	movq xmm1, rdi
	mulsd xmm0, xmm1
	mov rsp, rbp
	pop rbp
	ret

%macro atanTaylor64 1		; calculates the arctan of a 64-bit floating-point number
				; via a Taylor series
				; modifies the following registers: rdi, rsi, xmm0 - xmm6
	movsd xmm1, %1
	mulsd xmm1, xmm1	
	movq rdi, xmm1
	mov rsi, 0x8000000000000000
	or rdi, rsi
	movq xmm1, rdi		; -x*x
	mov rdi, 0x3ff0000000000000
	movq xmm3, rdi		; denominator
	mov rdi, 0x4000000000000000
	movq xmm4, rdi		; incrementor for denominator = 2.0
	movsd xmm5, %1		; power
	xorpd xmm6, xmm6	; test sum
				; %1 is the sum
	.taylorLoop:
		mulsd xmm5, xmm1
		movsd xmm2, xmm5
		addsd xmm3, xmm4
		divsd xmm2, xmm3
		addsd %1, xmm2
		ucomisd %1, xmm6
		movsd xmm6, %1
		jnz .taylorLoop
%endmacro

atan64:				; calculates the arctangent of a 64-bit floating-point number
				; modifies the following registers: rdi, rsi, rcx, rax, xmm0 - xmm6
	mov rcx, 0x8000000000000000
	movq rdi, xmm0
	and rcx, rdi
	xor rdi, rcx
	movq xmm0, rdi		; rcx has the information about the sign of the input
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi
	movsd xmm2, xmm0
	mulsd xmm2, xmm2
	addsd xmm2, xmm1
	sqrtsd xmm2, xmm2
	addsd xmm2, xmm1
	divsd xmm0, xmm2
	ucomisd xmm0, xmm1
	lahf
	jc .smallerThanOne
	
	divsd xmm1, xmm0
	movsd xmm0, xmm1
	
	.smallerThanOne:
		atanTaylor64 xmm0
		sahf
		jc .noRecip
		
		mov rdi, 0x3ff921fb54442d18
		movq xmm1, rdi
		subsd xmm1, xmm0
		movsd xmm0, xmm1
	
	.noRecip:
		movq rdi, xmm0
		xor rdi, rcx
		movq xmm0, rdi
		addsd xmm0, xmm0
	ret

asin64:				; calculates the arcsine of a 64-bit floating-point number
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	mov rdi, 0x3ff0000000000000
	movq xmm1, rdi
	movsd xmm2, xmm0
	mulsd xmm2, xmm2
	subsd xmm1, xmm2
	movq rdi, xmm1
	and rdi, rdi
	jz .one
	sqrtsd xmm1, xmm1
	divsd xmm0, xmm1
	
	call atan64
	
	mov rsp, rbp
        pop rbp
        ret
	
	.one:
		movq rdi, xmm0
                mov rsi, 0x3ff921fb54442d18
                mov rdx, 0xbff921fb54442d18
                mov rcx, 0x8000000000000000
                and rcx, rdi
                cmovnz rsi, rdx
	
	mov rsp, rbp
	pop rbp
	ret

acos64:				; calculates the arccosine of a 64-bit floating-point number
				; modifies the following registers: rdi, rsi, rax, rcx, rdx, xmm0 - xmm6
	push rbp
	mov rbp, rsp
	mov rdi, 0x3ff0000000000000
        movq xmm1, rdi
        movsd xmm2, xmm0
        mulsd xmm2, xmm2
        subsd xmm1, xmm2
	movq rdi, xmm1
	and rdi, rdi
	jz .one
        sqrtsd xmm1, xmm1
        divsd xmm0, xmm1

        call atan64
	
	mov rdi, 0x3ff921fb54442d18
	movq xmm1, rdi
	subsd xmm1, xmm0
	movsd xmm0, xmm1
	
	mov rsp, rbp
	pop rbp
	ret
	
	.one:
		movq rdi, xmm0
		xor rsi, rsi
		mov rdx, 0x400921fb54442d18
		mov rcx, 0x8000000000000000
		and rcx, rdi
		cmovnz rsi, rdx
		movq xmm0, rsi
		
	mov rsp, rbp
	pop rbp
	ret

agm64:				; calculates the arithmetic-geometric mean of two 64-bit floating-point numbers
				; modifies the following registers: rdi, rsi, xmm0 - xmm4
	mov rdi, 0x3fe0000000000000
	movq xmm4, rdi
	xorpd xmm3, xmm3
	
	.loop:
		movsd xmm2, xmm0
		addsd xmm0, xmm1
		mulsd xmm0, xmm4
		sqrtsd xmm2, xmm2	; do two sqrts to avoid over- or underflows due to the multiplication
		sqrtsd xmm1, xmm1
		mulsd xmm1, xmm2
		ucomisd xmm0, xmm3
		movsd xmm3, xmm0
		jnz .loop

	addsd xmm0, xmm1
	mulsd xmm0, xmm4
	
	ret

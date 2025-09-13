---
title: Coupled Cluster理论学习笔记（3）：振幅方程的推导（以 CCD 为例）
published: 2025-09-13
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---

&emsp;&emsp;与能量方程不同，振幅方程需要计算激发态矢量与有效哈密顿和真空态矢量之间的内积，所以振幅方程的推导要比能量方程复杂很多，但是具体推导过程类似，依旧是先把振幅方程自然截断再对剩下的非零项逐项进行收缩。。本文将以 CCD （即激发算符只保留双激发部分）为例，简要介绍 CC 振幅方程的推导过程。

&emsp;&emsp;在推导开始之前，这里先将哈密顿算符的各分子轨道指标换成占据/非占据轨道指标。这样做能够清晰地看到各项的 q-产生湮灭算符个数，以便后续推导。

$$
\begin{equation}
\begin{aligned}
f_{pq}\{a_p^\dagger a_q\}
&= f_{ij} \{a_i^\dagger a_j\} + f_{ab} \{a_a^\dagger a_b\} \\
&= f_{ia} \{a_i^\dagger a_a\} + f_{ai} \{a_a^\dagger a_i\} \\
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
\frac{1}{4} \braket{pq||rs} \{a_p^\dagger a_q^\dagger a_s a_r\}
&= \frac{1}{4} \braket{ab||cd} \{a_a^\dagger a_b^\dagger a_d a_c\} + \frac{1}{4} \braket{ij||kl} \{a_i^\dagger a_j^\dagger a_k a_l\} +  \braket{ia||bj} \{a_i^\dagger a_a^\dagger a_j a_b\} \\
&+ \frac{1}{2} \braket{ai||bc} \{a_a^\dagger a_i^\dagger a_c a_b\} + \frac{1}{2} \braket{ij||ka} \{a_i^\dagger a_j^\dagger a_k a_a\} + \frac{1}{2} \braket{ab||ci} \{a_a^\dagger a_b^\dagger a_i a_c\} \\
&+ \frac{1}{2} \braket{ia||jk} \{a_i^\dagger a_a^\dagger a_k a_j\} + \frac{1}{4} \braket{ij||ab} \{a_i^\dagger a_j^\dagger a_b a_a\} + \frac{1}{4} \braket{ab||ij} \{a_a^\dagger a_b^\dagger a_j a_i\}
\end{aligned}
\end{equation}
$$

其中两体算符部分使用了双电子积分的对称性，合并了一些项。可以看到两体算符前三项均包含两个产生算符两个湮灭算符、第四项和第五项三产生一湮灭、第六项和第七项一产生三湮灭、第八项四湮灭、最后一项四产生。

## 1. CCD 振幅方程的推导

&emsp;&emsp;在本系列第一篇博客我们已经推导过 CC 振幅方程。令激发算符 $\hat{T} = \hat{T}_2$ ，可以得到 CCD 双激发振幅所满足的方程，即

$$
\begin{equation}
\begin{aligned}
\Omega_{abij}
&= \bra{\Psi_{ij}^{ab}}e^{-\hat{T}}\hat{H} e^{\hat{T}}\ket{\Psi_0} \\
&= \bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}e^{-\hat{T}_2}\hat{H} e^{\hat{T}_2}\ket{\Psi_0} \\
&= 0
\end{aligned}
\end{equation}
$$

由于激发算符全为产生算符，而 $\{a^\dagger_i a^\dagger_j a_b a_a\}$ 为湮灭算符，为了保证最后的收缩表达式中能含有全收缩项，产生算符和湮灭算符的数量必需相等。由于哈密顿算符中的项最多包含四个湮灭算符，所以产生算符也最多只有八个。故上式最多只能保留两个双激发算符，因此上式的 BCH 展开式二阶以上的项全部恒为零，即

$$
\begin{equation}
\begin{aligned}
\Omega_{abij}
&= \bra{\Psi_{ij}^{ab}}e^{-\hat{T}_2}\hat{H} e^{\hat{T}_2}\ket{\Psi_0} \\
&= \bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}(\hat{H} + \hat{H}\hat{T}_2 + \frac{1}{2!}\hat{H}\hat{T}_2\hat{T}_2)\ket{\Psi_0} \\
&= 0
\end{aligned}
\end{equation}
$$

与能量方程的推导类似，接下来我们来展开计算上面这三项，以得到双激发振幅所满足的代数方程。下面推导过程若没有特殊说明，均使用爱因斯坦求和约定。

**（1）$\bra{\Psi_{ij}^{ab}}\hat{H}\ket{\Psi_0}$**

&emsp;&emsp;该项中含有四个湮灭算符，哈密顿算符中只有含有四个产生算符的项才能够全收缩，即

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}\hat{H}\ket{\Psi_0}
&=  \bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}\frac{1}{4} \braket{cd||kl} \{a_c^\dagger a_d^\dagger a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{4} (1-\hat{P}_{ab})(1-\hat{P}_{ij})\braket{ab||ij} \\
&= \braket{ab||ij}
\end{aligned}
\end{equation}
$$

其中 $\hat{P}_{ij}$ 为交换算符，其定义为 $\hat{P}_{ij} A_{abij} = A_{abji}$ 。

**（2）$\bra{\Psi_{ij}^{ab}}\hat{H}\hat{T}_2\ket{\Psi_0}$**

&emsp;&emsp;该项中含有四个产生算符四个湮灭算符，故哈密顿算符中只有产生算符和湮灭算符数量相等的项才能够全收缩，即

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}\hat{H}\hat{T}_2\ket{\Psi_0}
&=  \bra{\Psi_{ij}^{ab}}(f_{mn} \{a_m^\dagger a_n\} + f_{ef} \{a_e^\dagger a_f\} + \frac{1}{4} \braket{ef||gh} \{a_e^\dagger a_f^\dagger a_h a_g\} \\
&+ \frac{1}{4} \braket{mn||op} \{a_m^\dagger a_n^\dagger a_p a_o\} +  \braket{me||fn} \{a_m^\dagger a_e^\dagger a_n a_f\}) \hat{T}_2\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

下面我们来逐项计算上式

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}f_{mn} \{a_m^\dagger a_n\}  \hat{T}_2\ket{\Psi_0}
&= \frac{1}{4}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}f_{mn} \{a_m^\dagger a_n\}  t_{kl}^{cd}\{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{4}f_{mn} t_{kl}^{cd}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\} \{a_m^\dagger a_n\}  \{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= -f_{ii} t_{ij}^{ab}-f_{jj}t_{ij}^{ab}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}f_{ef} \{a_e^\dagger a_f\}  \hat{T}_2\ket{\Psi_0}
&= \frac{1}{4}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}f_{ef} \{a_e^\dagger a_f\} t_{kl}^{cd}\{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{4}f_{ef} t_{kl}^{cd}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\} \{a_e^\dagger a_f\}  \{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= f_{aa} t_{ij}^{ab} + f_{bb}t_{ij}^{ab}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}\frac{1}{4} \braket{ef||gh} \{a_e^\dagger a_f^\dagger a_h a_g\}  \hat{T}_2\ket{\Psi_0}
&= \frac{1}{4}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}\frac{1}{4} \braket{ef||gh} \{a_e^\dagger a_f^\dagger a_h a_g\} t_{kl}^{cd}\{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{16} \braket{ef||gh} t_{kl}^{cd}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\} \{a_e^\dagger a_f^\dagger a_h a_g\} \{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{16} (1-\hat{P}_{ij})(1-\hat{P}_{ab})(1-\hat{P}_{cd})\braket{ab||cd} t_{ij}^{cd} \\
&= \frac{1}{2} \braket{ab||cd} t_{ij}^{cd}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}\frac{1}{4} \braket{mn||op} \{a_m^\dagger a_n^\dagger a_p a_o\}   \hat{T}_2\ket{\Psi_0}
&= \frac{1}{4}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}\frac{1}{4} \braket{mn||op} \{a_m^\dagger a_n^\dagger a_p a_o\} t_{kl}^{cd}\{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{16} \braket{mn||op} t_{kl}^{cd}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\} \{a_m^\dagger a_n^\dagger a_p a_o\}  \{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{16} (1-\hat{P}_{ij})(1-\hat{P}_{ab})(1-\hat{P}_{kl})\braket{ij||kl} t_{kl}^{ab} \\
&= \frac{1}{2} \braket{ij||kl} t_{kl}^{ab}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}} \braket{me||fn} \{a_m^\dagger a_e^\dagger a_n a_f\} \hat{T}_2\ket{\Psi_0}
&= \frac{1}{4}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}\braket{me||fn} \{a_m^\dagger a_e^\dagger a_n a_f\} t_{kl}^{cd}\{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{4} \braket{me||fn} t_{kl}^{cd}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\} \{a_m^\dagger a_e^\dagger a_n a_f\}  \{a^\dagger_c a^\dagger_d a_l a_k\} \ket{\Psi_0} \\
&= \frac{1}{4} (1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kb||cj}[(1-\hat{P}_{ik})(1-\hat{P}_{ac}) t_{ik}^{ac}] \\
&=  (1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kb||cj} t_{ik}^{ac}
\end{aligned}
\end{equation}
$$

综上所述

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_{ij}^{ab}}\hat{H}\hat{T}_2\ket{\Psi_0}
&=  (f_{aa} + f_{bb} - f_{ii} - f_{jj} )t_{ij}^{ab} + \frac{1}{2} \braket{ab||cd} t_{ij}^{cd} \\
&+ \frac{1}{2}\braket{ij||kl} t_{kl}^{ab} + (1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kb||cj} t_{ik}^{ac}
\end{aligned}
\end{equation}
$$

**（3）$\frac{1}{2!}\bra{\Psi_{ij}^{ab}}\hat{H}\hat{T}_2\hat{T}_2\ket{\Psi_0}$**

&emsp;&emsp;该项包含四个湮灭算符和八个产生算符，所以只有哈密顿算符中包含四个湮灭算符的那一项才能够全收缩，即

$$
\begin{equation}
\begin{aligned}
\frac{1}{2!}\bra{\Psi_{ij}^{ab}}\hat{H}\hat{T}_2\hat{T}_2\ket{\Psi_0}
&=  \frac{1}{2!}\bra{\Psi_{ij}^{ab}}\frac{1}{4} \braket{op||gh} \{a_o^\dagger a_p^\dagger a_h a_g\}\hat{T}_2\hat{T}_2\ket{\Psi_0} \\
&= \frac{1}{128}\bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\} \braket{op||gh} \{a_o^\dagger a_p^\dagger a_h a_g\}t_{kl}^{cd}\{a^\dagger_c a^\dagger_d a_l a_d\}  t_{mn}^{ef}\{a^\dagger_e a^\dagger_f a_n a_n\} \ket{\Psi_0} \\
&= \frac{1}{128} \braket{op||gh} t_{kl}^{cd} t_{mn}^{ef} \bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}  \{a_o^\dagger a_p^\dagger a_h a_g\}\{a^\dagger_c a^\dagger_d a_l a_k\}  \{a^\dagger_e a^\dagger_f a_n a_m\} \ket{\Psi_0} \\
&= \frac{1}{128}\braket{kl||cd} ((1-\hat{P}_{ij})(1-\hat{P}_{ab})(1-\hat{P}_{cd})(1-\hat{P}_{kl})t_{ij}^{cd}t_{kl}^{ab} \\
&+ (1-\hat{P}_{ij})(1-\hat{P}_{ab})(1-\hat{P}_{cd})(1-\hat{P}_{kl})(1-\hat{P}_{ac})(1-\hat{P}_{bd})(1-\hat{P}_{ik})(1-\hat{P}_{jl})t_{ik}^{ac}t_{jl}^{bd} \\
&-  2(1-\hat{P}_{ij})(1-\hat{P}_{ab})(1-\hat{P}_{cd})(1-\hat{P}_{kl})(1-\hat{P}_{ik})(1-\hat{P}_{jl})t_{ik}^{ab}t_{jl}^{cd} \\
&- 2(1-\hat{P}_{ij})(1-\hat{P}_{ab})(1-\hat{P}_{cd})(1-\hat{P}_{kl})(1-\hat{P}_{ac})(1-\hat{P}_{bd})t_{ij}^{ac}t_{kl}^{bd})  \\
&= \frac{1}{4}\braket{kl||cd}t_{ij}^{cd}t_{kl}^{ab} + \frac{1}{2}(1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kl||cd}t_{ik}^{ac}t_{jl}^{bd} \\
&- \frac{1}{2}(1-\hat{P}_{ij})\braket{kl||cd}t_{ik}^{ab}t_{jl}^{cd} - \frac{1}{2}(1-\hat{P}_{ab})\braket{kl||cd}t_{ij}^{ac}t_{kl}^{bd}
\end{aligned}
\end{equation}
$$

## 2. 结论

&emsp;&emsp;综上所述，CCD 双激发振幅所满足的代数方程为

$$
\begin{equation}
\begin{aligned}
\Omega_{abij}
&= \bra{\Psi_{ij}^{ab}}e^{-\hat{T}_2}\hat{H} e^{\hat{T}_2}\ket{\Psi_0} \\
&= \bra{\Psi_0}\{a^\dagger_i a^\dagger_j a_b a_a\}(\hat{H} + \hat{H}\hat{T}_2 + \frac{1}{2!}\hat{H}\hat{T}_2\hat{T}_2)\ket{\Psi_0} \\
&= \braket{ab||ij} + (f_{aa} + f_{bb} - f_{ii} - f_{jj} )t_{ij}^{ab} + \frac{1}{2} \braket{ab||cd} t_{ij}^{cd} \\
&+ \frac{1}{2}\braket{ij||kl} t_{kl}^{ab} + (1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kb||cj} t_{ik}^{ac} \\
&+ \frac{1}{4}\braket{kl||cd}t_{ij}^{cd}t_{kl}^{ab} + \frac{1}{2}(1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kl||cd}t_{ik}^{ac}t_{jl}^{bd} \\
&- \frac{1}{2}(1-\hat{P}_{ij})\braket{kl||cd}t_{ik}^{ab}t_{jl}^{cd} - \frac{1}{2}(1-\hat{P}_{ab})\braket{kl||cd}t_{ij}^{ac}t_{kl}^{bd} \\
&= 0
\end{aligned}
\end{equation}
$$

通过求解上述方程组，即可得到双激发振幅 $t_{ij}^{ab}$。在实际计算中，迭代方程往往比代数方程更容易处理。将上式写成迭代方程的形式，令 $\varepsilon_p = f_{pp}$ 为分子轨道能量，可得

$$
\begin{equation}
\begin{aligned}
 t_{ij}^{ab}
&= \frac{1}{ \varepsilon_{i} + \varepsilon_{j} - \varepsilon_{a} - \varepsilon_{b} }\bigg[\braket{ab||ij}  + \frac{1}{2} \braket{ab||cd} t_{ij}^{cd} \\
&+ \frac{1}{2}\braket{ij||kl} t_{kl}^{ab} + (1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kb||cj} t_{ik}^{ac} \\
&+ \frac{1}{4}\braket{kl||cd}t_{ij}^{cd}t_{kl}^{ab} + \frac{1}{2}(1-\hat{P}_{ij})(1-\hat{P}_{ab})\braket{kl||cd}t_{ik}^{ac}t_{jl}^{bd} \\
&- \frac{1}{2}(1-\hat{P}_{ij})\braket{kl||cd}t_{ik}^{ab}t_{jl}^{cd} - \frac{1}{2}(1-\hat{P}_{ab})\braket{kl||cd}t_{ij}^{ac}t_{kl}^{bd} \bigg]
\end{aligned}
\end{equation}
$$

这就是一般量化代码中最常用的 CCD 方程的形式，可直接用来编写程序。若令方程右边的双激发振幅等于零，可得

$$
\begin{equation}
\begin{aligned}
 t_{ij}^{ab}
&= \frac{\braket{ab||ij} }{ \varepsilon_{i} + \varepsilon_{j} - \varepsilon_{a} - \varepsilon_{b} }
\end{aligned}
\end{equation}
$$

这就是 MP2 振幅，将其代入 CC 能量方程即可得到 MP2 能量。因此我们常使用 MP2 振幅作为双激发振幅的初猜来开始振幅方程的迭代。

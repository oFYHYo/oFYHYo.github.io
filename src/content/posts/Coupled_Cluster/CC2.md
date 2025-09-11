---
title: Coupled Cluster理论学习笔记（2）：能量方程的推导
published: 2025-09-11
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---

&emsp;&emsp;半年前开的坑到现在我还没开始填。趁着现在没啥事，多写几篇 CC 教程。

&emsp;&emsp;本文主要介绍 Coupled Cluster 的能量方程的推导。根据相关能的定义以及上篇博客的结论，可以得到 CC 能量的表达式，即

$$
\begin{equation}
\begin{aligned}
E &= \bra{\Psi_0}e^{-\hat{T}}\hat{H}e^{\hat{T}}\ket{\Psi_0} \\
&= \bra{\Psi_0}\hat{H}_{cc}\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

其中 $\hat{H}_{cc}=e^{-\hat{T}}\hat{H}e^{\hat{T}}$ 为 CC 的有效哈密顿量，$\hat{T} = \sum\hat{T}_n$ 为激发算符。各符号的详细定义请参考上一篇博客。下面我将会从这个公式出发，推导出实际用于计算的，用振幅表示的能量表达式。

## 1. 哈密顿算符的正则序形式

&emsp;&emsp;为了简化后续推导，这里先导出哈密顿算符的正则序形式。由于正则序一般为左产生右湮灭的形式，故正则序在真空态下的期望为0，所以将哈密顿算符写成正则序形式可以大大简化后续的推导过程。本章将会从二次量子化形式的哈密顿量出发，导出哈密顿的正则序形式，即

$$
\begin{equation}
\hat{H} = E_{HF} + \sum_{pq} f_{pq} \{a_p^\dagger a_q\} + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} \{a_p^\dagger a_q^\dagger a_s a_r\}
\end{equation}
$$

其中 $f_{pq} = h_{pq} +\sum_i\bra{pi}\ket{qi}$ 为体系的 fock 矩阵，$E_{HF} = \sum_{i} h_{ii} +\frac{1}{2}\sum_{ij} \bra{ij}\ket{ij}$ 为 HF 基态（参考态）能量。

&emsp;&emsp;不考虑相对论，在BO近似下，分子体系的哈密顿算符由两部分组成，单体算符（单电子哈密顿算符$\hat{h}$）和两体算符（电子排斥算符$\frac{1}{r_{12}}$），可参考树桑的文章[^1]写出哈密顿算符的二次量子化形式，即

$$
\begin{equation}
\begin{aligned}
\hat{H} &= \sum_{pq} h_{pq} a_p^\dagger a_q + \frac{1}{2}\sum_{pqrs} \braket{pq|rs} a_p^\dagger a_q^\dagger a_s a_r \\
&= \sum_{pq} h_{pq} a_p^\dagger a_q + \frac{1}{4}\sum_{pqrs} \braket{pq|rs} a_p^\dagger a_q^\dagger a_s a_r + \frac{1}{4}\sum_{pqrs} \braket{pq|sr} a_p^\dagger a_q^\dagger a_r a_s \\
&= \sum_{pq} h_{pq} a_p^\dagger a_q + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} a_p^\dagger a_q^\dagger a_s a_r
\end{aligned}
\end{equation}
$$

其中第二个等号运用了费米子湮灭算符的交换反对称性，即 $a_r a_s = - a_s a_r$。运用 wick 定理，可将哈密顿算符用正则序表示出来，详细推导过程如下

$$
\begin{equation}
\begin{aligned}
\hat{H} &= \sum_{pq} h_{pq} a_p^\dagger a_q + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} a_p^\dagger a_q^\dagger a_s a_r \\
&= \sum_{pq} h_{pq} \{a_p^\dagger a_q\} + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} \{a_p^\dagger a_q^\dagger a_s a_r\} \\
&+\sum_{pq} h_{pq} \{a_p^{\dagger\bullet} a_q^\bullet\} \\
 &+\frac{1}{4}\sum_{pqrs} \braket{pq||rs} (\{a_p^{\dagger\bullet} a_q^\dagger a_s^\bullet a_r\}+\{a_p^\dagger a_q^{\dagger\bullet} a_s^\bullet a_r\}+\{a_p^{\dagger\bullet} a_q^\dagger a_s a_r^\bullet\}+\{a_p^\dagger a_q^{\dagger\bullet} a_s a_r^\bullet\}) \\
 &+\frac{1}{4}\sum_{pqrs} \braket{pq||rs}(\{a_p^{\dagger\bullet} a_q^{\dagger\bullet\bullet} a_s^\bullet a_r^{\bullet\bullet}\}+\{a_p^{\dagger\bullet\bullet} a_q^{\dagger\bullet} a_s^\bullet a_r^{\bullet\bullet}\}) \\
&= \sum_{pq} h_{pq} \{a_p^\dagger a_q\} + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} \{a_p^\dagger a_q^\dagger a_s a_r\} \\
&+\sum_{pq} h_{pq} \delta_{pq}\delta_{pi} +\frac{1}{4}\sum_{pqrs} \braket{pq||rs}(\delta_{pr}\delta_{pi}\delta_{qs}\delta_{qj} - \delta_{ps}\delta_{pi}\delta_{qr}\delta_{qj})\\
&+\frac{1}{4}\sum_{pqrs} \braket{pq||rs} (-\{a_q^\dagger a_r\}\delta_{ps}\delta_{pi}+\{a_p^\dagger a_r\}\delta_{qs}\delta_{qi}+\{a_q^\dagger a_s \}\delta_{pr}\delta_{pi}-\{a_p^\dagger a_s \}\delta_{qr}\delta_{qi}) \\
&= \sum_{pq} (h_{pq} +\sum_i\braket{pi||qi})\{a_p^\dagger a_q\} + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} \{a_p^\dagger a_q^\dagger a_s a_r\}  \\
&+\sum_{i} h_{ii} +\frac{1}{2}\sum_{ij} \braket{ij||ij}\\
&= E_{HF} + \sum_{pq} f_{pq} \{a_p^\dagger a_q\} + \frac{1}{4}\sum_{pqrs} \braket{pq||rs} \{a_p^\dagger a_q^\dagger a_s a_r\}
\end{aligned}
\end{equation}
$$

这就是分子体系哈密顿算符的正则序形式。由于正则序为左产生右湮灭的形式，与零阶波函数（即费米真空态）的内积为零，容易验证，$\bra{\Psi_0}\hat{H} \ket{\Psi_0}=E_{HF}$，这也是 HF 基态能量的最初定义。接下来我们将使用哈密顿量的正则序形式进行后续的推导。

## 2. 能量展开式的自然截断

&emsp;&emsp;将 CC 有效哈密顿用上篇文章提到的 BCH 公式展开，可得

$$
\begin{equation}
\begin{aligned}
\hat{H}_{cc}&=e^{-\hat{T}}\hat{H}e^{\hat{T}}\\
&= \hat{H}+[\hat{H},\hat{T}]+\frac{1}{2!}[[\hat{H},\hat{T}],\hat{T}]+\frac{1}{3!}[[[\hat{H},\hat{T}],\hat{T}],\hat{T}]+\cdots
\end{aligned}
\end{equation}
$$

看上去要计算无穷项，实则不然。在试剂计算能量的时候，通常要将有效哈密顿在 HF 基态上求期望，而正则序部分的期望为零，所以只需要找到上式的非正则序部分并计算期望就能得到最终的能量表达式。实际上，由于自然截断，上式只有四项存在非正则序部分，即

$$
\begin{equation}
\begin{aligned}
E_{CC} &= E_{HF} + \bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0} + \bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0} + \frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

第一项很好理解，因为 HF 基态能量是常数项。接下来我们将会讨论为什么只有后面几项存在非正则序部分。

&emsp;&emsp;下面的推导过程需要计算正则序之间的乘积。正则序之间的乘积同样可以使用wick定理进行计算，这里有一点需要注意。正则序是一种左产生右湮灭的算符形式，本质上和普通算符一样，可以表示为一系列产生湮灭算符的乘积。由于正则序内部的收缩只可能是产生产生、湮灭湮灭和左产生右湮灭这三种情况，而这三种收缩都为零，所以正则序内部的收缩都为零。所以在应用 wick 定理计算正则序之间的乘积时，只有正则序之间的收缩非零。

&emsp;&emsp;在第一部分的推导中我们可以知道，哈密顿算符由单体算符、两体算符和常数项组成，为了使得哈密顿与激发算符的乘积不为零，根据 wick 定理和正则序的定义，只有当哈密顿算符与激发算符的乘积存在完全收缩，才有可能得到非正则序项，例如

$$
\begin{equation}
\begin{aligned}
\{a_p^\dagger a_q\}\{a_a^\dagger a_i\} \rightarrow \{a_p^{\dagger\bullet} a_q^{\bullet\bullet}\}\{a_a^{\dagger\bullet\bullet} a_i^\bullet\} = \delta_{pi}\delta_{qa}
\end{aligned}
\end{equation}
$$

上述正则序的乘积中所有的产生湮灭算符都得到了收缩，最后得到了不含正则序的项。而下面的这些乘积则不能完全收缩，例如

$$
\begin{equation}
\begin{aligned}
\{a_p^\dagger a_q\}\{a_a^\dagger a_i\} \{a_b^\dagger a_j\} , \{a_p^\dagger a^\dagger_q a_s a_r\}\{a_a^\dagger a^\dagger_b a_j a_i\} \{a_c^\dagger a_k\} \cdots
\end{aligned}
\end{equation}
$$

这些项最后收缩得到的结果都含有正则序，所以这些项作用到真空态上会等于零。由于哈密顿算符只包含单体部分和两体部分，所以所有包含 $\hat{T}_3$ 及以上激发算符的项都为零；二阶以上展开项都为零，即

$$
\begin{equation}
\begin{aligned}
E_{CC} &= E_{HF} + \bra{\Psi_0}[\hat{H}, \hat{T}_1+\hat{T}_2]\ket{\Psi_0} + \frac{1}{2!}\bra{\Psi_0}[[\hat{H}, \hat{T}_1] ,\hat{T_1}]\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

又因为激发算符都是 q-产生算符，而 q-产生算符在左边的收缩表达式为零，所以只有激发算符在哈密顿算符右边的项非零。简化上述对易表达式，最后得到

$$
\begin{equation}
\begin{aligned}
E_{CC} &= E_{HF} + \bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0} + \bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0} + \frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

这就是经过自然截断保留非零项的 CC 能量表达式。第一项为 HF 基态能量，后面三项为 **CC 相关能** 。现在我们将原来的能量表达式的无穷求和给简化成上述四项。在下一部分我们将详细计算上述几项，得到能用振幅表示的能量表达式。

## 3. 能量方程的推导

&emsp;&emsp;上一部分中我们得到了自然截断后的能量表达式，接下来我们逐项来进行计算。为了简化推导过程，下方推导均使用了爱因斯坦求和约定。

**（1） $\bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0}$**

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0} &= \bra{\Psi_0}f_{pq} t_i^a\{a_p^\dagger a_q\}\{a^\dagger_a a_i\}\ket{\Psi_0} \\
&= \bra{\Psi_0}f_{pq} t_i^a\{a_p^{\dagger\bullet\bullet} a_q^\bullet\}\{a^{\dagger\bullet}_a a_i^{\bullet\bullet}\}\ket{\Psi_0} \\
&= \bra{\Psi_0}f_{pq} t_i^a \delta_{pi} \delta_{qa}\ket{\Psi_0} \\
&= f_{ia} t_i^a
\end{aligned}
\end{equation}
$$

**（2） $\bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0}$**

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0} &= \frac{1}{16}\bra{\Psi_0} \braket{pq||rs} t_{ij}^{ab}\{a^\dagger_p a^\dagger_q a_s a_r\} \{a^\dagger_a a^\dagger_b a_j a_i\}\ket{\Psi_0} \\
&= \frac{1}{16} \braket{pq||rs} t_{ij}^{ab}(\delta_{pi}\delta_{qj}\delta_{ra}\delta_{sb} - \delta_{pi}\delta_{qj}\delta_{rb}\delta_{sa} - \delta_{pj}\delta_{qi}\delta_{ra}\delta_{sb} + \delta_{pj}\delta_{qi}\delta_{rb}\delta_{sa})\\
&= \frac{1}{16}  t_{ij}^{ab}(\braket{ij||ab}- \braket{ij||ba} - \braket{ji||ab} + \braket{ji||ba})\\
&= \frac{1}{4}  t_{ij}^{ab}\braket{ij||ab}
\end{aligned}
\end{equation}
$$

**（3） $\frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0}$**

$$
\begin{equation}
\begin{aligned}
\frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0} &= \frac{1}{8}\bra{\Psi_0} \braket{pq||rs} t_i^a t_j^b\{a^\dagger_p a^\dagger_q a_s a_r\} \{a^\dagger_a a_i \}\{a^\dagger_b a_j \} \ket{\Psi_0} \\
&= \frac{1}{8}\braket{pq||rs} t_i^a t_j^b(\delta_{pi}\delta_{qj}\delta_{ra}\delta_{sb} - \delta_{pi}\delta_{qj}\delta_{rb}\delta_{sa} - \delta_{pj}\delta_{qi}\delta_{ra}\delta_{sb} + \delta_{pj}\delta_{qi}\delta_{rb}\delta_{sa}) \\
&= \frac{1}{8} t_i^a t_j^b(\braket{ij||ab}- \braket{ij||ba} - \braket{ji||ab} + \braket{ji||ba}) \\
&= \frac{1}{2}  t_i^a t_j^b\braket{ij||ab}
\end{aligned}
\end{equation}
$$

&emsp;&emsp;综上所述，CC 相关能表达式为

$$
\begin{equation}
\begin{aligned}
E_{corr} &= \bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0} + \bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0} + \frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0} \\
&= f_{ia} t_i^a + \frac{1}{4}  t_{ij}^{ab}\braket{ij||ab} + \frac{1}{2}  t_i^a t_j^b\braket{ij||ab}
\end{aligned}
\end{equation}
$$

将求和号写出，加上 HF 基态能量，即可得到耦合簇理论的能量表达式，即

$$
\begin{equation}
E_{CC} = E_{HF} + \sum_{ia}f_{ia} t_i^a + \frac{1}{4} \sum_{ijab} \braket{ij||ab}(t_{ij}^{ab} + 2 t_i^a t_j^b)
\end{equation}
$$

可以看到，由于自然截断的存在，无论激发算符展开到多少阶，能量表达式中始终只含有单激发振幅和双激发振幅，并不直接包含更高阶的激发振幅。所以一般来说 CC 展开到双激发（CCSD）就已经足够精确了，展开到更高阶的激发算符也只是对单双激发振幅加上修正而已。

[^1]: <https://zhuanlan.zhihu.com/p/97685686>

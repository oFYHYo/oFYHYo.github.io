---
title: Coupled Cluster理论学习笔记（2）：能量方程的推导
published: 2025-09-10
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---

&emsp;&emsp;半年前开的坑到现在我还没开始填。趁着现在没啥事，多写几篇 CC 教程。

&emsp;&emsp;本文主要介绍 Coupled Cluster 的能量方程的推导。根据相关能的定义以及上篇博客的结论，可以得到相关能的表达式，即

$$
\begin{equation}
E_{corr} = \bra{\Psi_0}\hat{H}_{cc} - E_{HF}\ket{\Psi_0}
\end{equation}
$$

其中 $\hat{H}_{cc}=e^{-\hat{T}}\hat{H}e^{\hat{T}}$ 为 CC 的有效哈密顿量，$\hat{T} = \sum\hat{T}_n$ 为激发算符。各符号的详细定义请参考上一篇博客。

## 1. 哈密顿算符的正则序形式

&emsp;&emsp;由于正则序一般为左产生右湮灭的形式，故正则序在真空态下的期望为0，所以将哈密顿算符写成正则序形式可以大大简化后续的推导过程。本章将会从二次量子化形式的哈密顿量出发，导出哈密顿的正则序形式，即

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

$$
\begin{equation}
\begin{aligned}
\hat{H}_{cc}&=e^{-\hat{T}}\hat{H}e^{\hat{T}}\\
&= \hat{H}+[\hat{H},\hat{T}]+\frac{1}{2!}[[\hat{H},\hat{T}],\hat{T}]+\frac{1}{3!}[[[\hat{H},\hat{T}],\hat{T}],\hat{T}]+\cdots
\end{aligned}
\end{equation}
$$

## 3. 能量方程的推导

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
&= \frac{1}{16}  t_{ij}^{ab}(\bra{ij}\ket{ab}- \bra{ij}\ket{ba} - \bra{ji}\ket{ab} + \bra{ji}\ket{ba})\\
&= \frac{1}{4}  t_{ij}^{ab}\bra{ij}\ket{ab}
\end{aligned}
\end{equation}
$$

**（3） $\frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0}$**

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0} &= \frac{1}{8}\bra{\Psi_0} \braket{pq||rs} t_i^a t_j^b\{a^\dagger_p a^\dagger_q a_s a_r\} \{a^\dagger_a a_i \}\{a^\dagger_b a_j \} \ket{\Psi_0} \\
&= \frac{1}{8}\braket{pq||rs} t_i^a t_j^b(\delta_{pi}\delta_{qj}\delta_{ra}\delta_{sb} - \delta_{pi}\delta_{qj}\delta_{rb}\delta_{sa} - \delta_{pj}\delta_{qi}\delta_{ra}\delta_{sb} + \delta_{pj}\delta_{qi}\delta_{rb}\delta_{sa}) \\
&= \frac{1}{8} t_i^a t_j^b(\bra{ij}\ket{ab}- \bra{ij}\ket{ba} - \bra{ji}\ket{ab} + \bra{ji}\ket{ba}) \\
&= \frac{1}{2}  t_i^a t_j^b\bra{ij}\ket{ab}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
E_{corr} &= \bra{\Psi_0}\hat{H}_{cc}\ket{\Psi_0} \\
&= \bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0} + \bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0} + \frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0} \\
&= f_{ia} t_i^a + \frac{1}{4}  t_{ij}^{ab}\bra{ij}\ket{ab} + \frac{1}{2}  t_i^a t_j^b\bra{ij}\ket{ab}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
E_{corr} = \sum_{ia}f_{ia} t_i^a + \frac{1}{4} \sum_{ijab} \bra{ij}\ket{ab}(t_{ij}^{ab} + 2 t_i^a t_j^b)
\end{equation}
$$

[^1]: <https://zhuanlan.zhihu.com/p/97685686>

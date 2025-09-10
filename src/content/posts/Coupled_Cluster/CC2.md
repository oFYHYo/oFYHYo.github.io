---
title: Coupled Cluster理论学习笔记（2）：能量方程的推导
published: 2025-09-091
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---

&emsp;&emsp;半年前开的坑到现在我还没开始填。正好最近没啥事，我来把CC教程给补充完。

&emsp;&emsp;本文主要介绍 Coupled Cluster 的能量方程的推导。

$$
\begin{equation}
E_{corr} = \bra{\Psi_0}\hat{H}_{cc}\ket{\Psi_0}
\end{equation}
$$

$$
\begin{equation}
\hat{H}_{cc}=e^{-\hat{T}}\hat{H}e^{\hat{T}}
\end{equation}
$$

### $\bra{\Psi_0}\hat{H} \hat{T}_1\ket{\Psi_0}$

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

### $\bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0}$

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_0}\hat{H} \hat{T}_2\ket{\Psi_0} &= \frac{1}{16}\bra{\Psi_0} \bra{pq}\ket{rs} t_{ij}^{ab}\{a^\dagger_p a^\dagger_q a_s a_r\} \{a^\dagger_a a^\dagger_b a_j a_i\}\ket{\Psi_0} \\
&= \frac{1}{16} \bra{pq}\ket{rs} t_{ij}^{ab}(\delta_{pi}\delta_{qj}\delta_{ra}\delta_{sb} - \delta_{pi}\delta_{qj}\delta_{rb}\delta_{sa} - \delta_{pj}\delta_{qi}\delta_{ra}\delta_{sb} + \delta_{pj}\delta_{qi}\delta_{rb}\delta_{sa})\\
&= \frac{1}{16}  t_{ij}^{ab}(\bra{ij}\ket{ab}- \bra{ij}\ket{ba} - \bra{ji}\ket{ab} + \bra{ji}\ket{ba})\\
&= \frac{1}{4}  t_{ij}^{ab}\bra{ij}\ket{ab}
\end{aligned}
\end{equation}
$$

### $\frac{1}{2!}\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0}$

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_0}\hat{H} \hat{T}_1 \hat{T_1}\ket{\Psi_0} &= \frac{1}{8}\bra{\Psi_0} \bra{pq}\ket{rs} t_i^a t_j^b\{a^\dagger_p a^\dagger_q a_s a_r\} \{a^\dagger_a a_i \}\{a^\dagger_b a_j \} \ket{\Psi_0} \\
&= \frac{1}{8}\bra{pq}\ket{rs} t_i^a t_j^b(\delta_{pi}\delta_{qj}\delta_{ra}\delta_{sb} - \delta_{pi}\delta_{qj}\delta_{rb}\delta_{sa} - \delta_{pj}\delta_{qi}\delta_{ra}\delta_{sb} + \delta_{pj}\delta_{qi}\delta_{rb}\delta_{sa}) \\
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

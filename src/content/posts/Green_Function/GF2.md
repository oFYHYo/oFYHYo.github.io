---
title: 格林函数学习笔记（2）：海森堡表象与场算符的运方程
published: 2025-12-05
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---


&emsp;&emsp;本文将主要介绍量子力学中的海森堡表象(Heisenberg picture，也可以叫海森堡绘景)以及海森堡表象下场算符的运动方程。

## 时间演化算符

&emsp;&emsp;态矢量$\ket{\Psi(t)}$所满足的薛定谔方程为

$$
\begin{equation}
    i\frac{d}{dt} \ket{\Psi(t)} = \hat{H} \ket{\Psi(t)}
\end{equation}
$$

&emsp;&emsp;定义时间演化算符为$\ket{\Psi(t)} = U(t,t_0)\ket{\Psi(t_0)},(t>t_0)$，代入薛定谔方程，可得时间演化算符所满足的方程为

$$
\begin{equation}
    i\frac{d}{dt} U(t,t_0) = \hat{H} U(t,t_0)
\end{equation}
$$

&emsp;&emsp;当$t=t_0$时，$\ket{\Psi(t)} = \hat{U}(t,t_0)\ket{\Psi(t_0)} = \ket{\Psi(t_0)}$故时间演化算符满足初值条件$\hat{U}(t_0,t_0) = \hat{\mathbf{1}}$。由于量子力学体系一般具有时间反演不变性，故时间演化算符是酉算符，定义$U(t_0,t) = U^\dagger(t,t_0)= U^{-1}(t,t_0)$，可得

$$
\begin{equation}
    -i\frac{d}{dt} U(t_0,t) = U(t_0,t) \hat{H}
\end{equation}
$$

&emsp;&emsp;对方程(1)两边积分，代入初值条件，可得

$$
\begin{equation}
    \hat{U}(t,t_0) = \mathbf{1} - i \int_{t_0}^{t} d t_1 \hat{H}(t_1) \hat{U}(t_1,t_0)
\end{equation}
$$

&emsp;&emsp;将演化算符$\hat{U}(t_1,t_0)$继续展开，可得演化算符的迭代级数

$$
\begin{equation}
    \begin{aligned}
        U(t,t_0)
        &= \mathbf{1} - i \int_{t_0}^{t} d t_1 \hat{H}(t_1) + (-i)^2 \int_{t_0}^{t} dt_1 \int_{t_0}^{t_1} dt_2 \hat{H}(t_1) \hat{H}(t_2) + \cdots \\
        &= \mathbf{1} + \sum_{n=1}^{\infty} (-i)^n \int_{t_0}^{t} dt_1 \int_{t_0}^{t_1} dt_2 \cdots \int_{t_0}^{t_{n-1}} dt_n \hat{H}(t_1) \hat{H}(t_2) \cdots \hat{H}(t_n) \\
        &= \mathbf{1} + \sum_{n=1}^{\infty} \frac{(-i)^n}{n!} \int_{t_0}^{t} dt_1 \int_{t_0}^{t} dt_2 \cdots \int_{t_0}^{t} dt_n T[\hat{H}(t_1) \hat{H}(t_2) \cdots \hat{H}(t_n)] \\
        &= T\bigg[e^{-i\int_{t_0}^{t} dz \hat{H}(z)}\bigg]
    \end{aligned}
\end{equation}
$$

&emsp;&emsp;其中$T$为时序算符，其定义为将算符按时间从大到小排列。第三个等号中的$n!$是因为在前面的积分中，时间变量$t_1,t_2,\cdots,t_n$的积分区间是有序的，而在最后一个等号中，时间变量$t_1,t_2,\cdots,t_n$的积分区间是无序的，故需要除以$n!$以消除多重计数。

&emsp;&emsp;故时间演化算符的一个形式解可以写成

$$
\begin{equation}
    U(t,t_0) = T[e^{-i\int_{t_0}^{t} dz \hat{H}(z)}]
\end{equation}
$$

&emsp;&emsp;若不同时间的哈密顿量对易，则上式中的时序算符$T$可以省略，$ U(t,t_0) = e^{-i\int_{t_0}^{t} dz \hat{H}(z)}$。若哈密顿量不显含时间，则时间演化算符可以进一步简化为$ U(t,t_0) = e^{-i \hat{H}(t-t_0)}$。

&emsp;&emsp;同理可得

$$
\begin{equation}
    U(t_0,t) = \bar{T}[e^{i\int_{t_0}^{t} dz \hat{H}(z)}]
\end{equation}
$$

&emsp;&emsp;其中$\bar{T}$为逆时序算符，其定义为将算符按时间从小到大排列。故可以得到时序算符的表达式

$$
\begin{equation}
    \hat{U}(t_1,t_2) =
    \begin{cases}
        T[e^{-i\int_{t_2}^{t_1} dz \hat{H}(z)}], & t_1 > t_2 \\
        \bar{T}[e^{i\int_{t_1}^{t_2} dz \hat{H}(z)}], & t_1 < t_2 \\
        \mathbf{1}, & t_1 = t_2
    \end{cases}
\end{equation}
$$

## 海森堡表象下场算符的运动方程

&emsp;&emsp;算符$\hat{O}(t)$在态矢量$\ket{\Psi(t)}$上的期望可以写成如下形式

$$
\begin{equation}
    \begin{aligned}
        O(t)
        &=\bra{\Psi(t)}\hat{O}(t)\ket{\Psi(t )}\\
        &=\bra{\Psi} U(t_0,t) \hat{O}(t) U(t,t_0)\ket{\Psi} = \bra{\Psi}  \hat{O}_H(t)\ket{\Psi}
    \end{aligned}
\end{equation}
$$

&emsp;&emsp;其中$\hat{O}_H(t) = U(t_0,t) \hat{O}(t) U(t,t_0) $为海森堡表象下算符的表达式。上式第一个等号代表薛定谔表象(schrödinger picture)下的算符与态矢量，第二个等号后代表海森堡表象(Heisenberg picture)下的算符与态矢量。若算符$\hat{O}$不显含时，可以看到，薛定谔表象下算符不随时间变化，态矢量随时间变化；而海森堡表象下算符随时间变化，态矢量不随时间变化。可以从时间演化算符的运动方程得到海森堡表象下算符的运动方程：

$$
\begin{equation}
    \begin{aligned}
        i\frac{d}{dt} \hat{O}_H(t)
        &= i\frac{d}{dt} \bigg(U(t_0,t) \hat{O}(t) U(t,t_0)\bigg) \\
        &= i\frac{d}{dt} U(t_0,t) \hat{O}(t) U(t,t_0) + U(t_0,t) i\frac{d}{dt} \hat{O}(t) U(t,t_0) + U(t_0,t) \hat{O}(t) i\frac{d}{dt} U(t,t_0) \\
        &= -U(t_0,t) \hat{H} \hat{O}(t) U(t,t_0) + U(t_0,t) i\frac{d}{dt} \hat{O}(t) U(t,t_0) + U(t_0,t) \hat{O}(t) \hat{H} U(t,t_0) \\
        &= U(t_0,t) \bigg(i\frac{d}{dt} \hat{O}(t) + [\hat{O}(t), \hat{H}] \bigg) U(t,t_0)  \\
        &= \bigg[\hat{O}_H(t),\hat{H}_H(t)\bigg] + i\frac{\partial}{\partial t} \hat{O}_H(t)
    \end{aligned}
\end{equation}
$$

&emsp;&emsp;其中$\frac{\partial}{\partial t} \hat{O}_H(t)= U(t_0,t)\frac{d}{dt} \hat{O}(t) U(t,t_0) $代表只对算符$\hat{O}(t)$求导而不对演化算符求导。若算符不显含时间，则上式简化为

$$
\begin{equation}
    i\frac{d}{dt} \hat{O}_H(t) = [\hat{O}_H(t),\hat{H}_H(t)]
\end{equation}
$$

&emsp;&emsp;可以看到，海森堡表象下算符的运动方程与经典力学中的泊松括号形式的运动方程形式上非常类似。如果算符与哈密顿量对易，则该算符在海森堡表象下不随时间变化，该算符对应的物理量就是一个守恒量。

&emsp;&emsp;在上一篇文章中我们推导了哈密顿量的二次量子化形式

$$
\begin{equation}
    \hat{H} = \int dx \, \psi^\dagger(x) \hat{h}(x)\psi(x) + \frac{1}{2} \int dx dx' \, \psi^\dagger(x) \psi^\dagger(x') v(x,x') \psi(x') \psi(x)
\end{equation}
$$

&emsp;&emsp;场算符在海森堡表象下的表达式为

$$
\begin{equation}
    \begin{gathered}
        \psi_{H}(x,t) = U(t_0,t) \psi(x) U(t,t_0) \\
        \psi_{H}^\dagger(x,t) = U(t_0,t) \psi^\dagger(x) U(t,t_0)
    \end{gathered}
\end{equation}
$$

&emsp;&emsp;将上式代入算符的运动方程(10)中，注意到场算符本身不显含时间，可以得到

$$
\begin{equation}
    \begin{aligned}
        i\frac{d}{dt} \psi_{H}(x,t)
        &= [\psi_{H}(t), \hat{H}_H] \\
        &= \int dx' \hat{h}[x'](\psi_{H}(x,t), \psi_{H}^\dagger(x') \psi_{H}(x')) \\
        &+ \frac{1}{2} \int dx' dx'' v(x',x'') [\psi_{H}(x,t), \psi_{H}^\dagger(x') \psi_{H}^\dagger(x'') \psi_{H}(x'') \psi_{H}(x')] \\
        &= \hat{h}\psi_{H}(x,t) + \int dx' v(x,x') \psi_{H}^\dagger(x',t) \psi_{H}(x',t) \psi_{H}(x,t)
    \end{aligned}
\end{equation}
$$

&emsp;&emsp;同理可得

$$
\begin{equation}
    -i\frac{d}{dt} \psi_{H}^\dagger(x,t)
    = \psi_{H}^\dagger(x,t) \hat{h} + \int dx' v(x,x') \psi_{H}^\dagger(x,t) \psi_{H}^\dagger(x',t) \psi_{H}(x',t)
\end{equation}
$$

式(14)和式(15)就是场算符的运动方程。

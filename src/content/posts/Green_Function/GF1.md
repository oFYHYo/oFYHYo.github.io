---
title: 格林函数学习笔记（1）：场算符与二次量子化
published: 2025-11-01
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---

## 前言

&emsp;&emsp;最近在学格林函数，记录一些学习笔记以备后续查阅。

&emsp;&emsp;本系列文章将主要讲述有限温度下多粒子格林函数（1）的性质、运动方程及其应用。

$$
\begin{equation}
    G(1,\cdots,n;1',\cdots,n')
    = \frac{1}{i^n}\frac{Tr\bigg[\mathcal{T}\big(e^{-i\int_{\gamma}dz \hat{H}(z)} \psi(1)\cdots\psi(n)\psi^\dagger(n')\cdots\psi^\dagger(1')\big)\bigg]}
    {Tr\bigg[\mathcal{T}\big(e^{-i\int_{\gamma}dz \hat{H}(z)}\big)\bigg]}
\end{equation}
$$

&emsp;&emsp;本文将会引入二次量子化，介绍场算符的概念并且用场算符表示哈密顿量。我的恩师在课上曾经说过：对一个体系量子化的步骤主要有以下两步：

1. 找到体系所在的希尔伯特空间
2. 找到体系所对应的希尔伯特空间中的哈密顿量

下文中将会主要通过这两步来引入二次量子化方法。本系列文章主要讨论费米子体系，因此后面若无特殊说明，所有讨论均假设粒子为费米子。

## 多粒子体系的量子力学

### 多粒子体系的希尔伯特空间

&emsp;&emsp;引入Fock空间$\mathcal{F}$，其定义为

$$
\begin{equation}
    \mathcal{F} = \{\mathcal{H}_0,\mathcal{H}_1,\cdots,\mathcal{H}_N,\cdots\}
\end{equation}
$$

其中$\mathcal{H}_N$为N粒子希尔伯特空间

&emsp;&emsp;令$\ket{\Phi_N} \in \mathcal{H}_N, \ket{\Phi_M} \in \mathcal{H}_M$其中$N \neq M$，定义不同粒子数的希尔伯特空间之间的内积为

$$
\begin{equation}
    \braket{\Phi_N}{\Phi_M} = 0
\end{equation}
$$

即Fock空间中源自不同粒子数的希尔伯特空间的态矢量要求正交。只有当两个矢量来自同一个多粒子希尔伯特空间时，上述内积才可能不为零。可以证明，上面定义的Fock空间满足希尔伯特空间的定义。之后若无特殊说明，默认体系所在的希尔伯特空间为Fock空间。

&emsp;&emsp;下面我们来定义一个特殊的态矢量——真空态。定义真空态为$\ket{0} \in \mathcal{H}_0$。由于$\mathcal{H}_0$中只含有这一个矢量，因此$\mathcal{H}_{0}$中的正交归一关系可以直接写成

$$
\begin{equation}
    \braket{0}{0} = 1
\end{equation}
$$

从上式可以知道，真空态矢量跟除了它自己以外的任意态矢量之间的内积为零。需要注意的是，真空态并不是Fock空间中的零矢量$\ket{\emptyset}$。真空态是一个有物理意义的态，而零矢量是非物理的。

### 多粒子希尔伯特空间中的基矢的一些基本性质

&emsp;&emsp;引入N粒子希尔伯特空间中的坐标基矢为$\ket{x_1\cdots x_N} \in \mathcal{H}_N$，其中$\{x_i\} = \{r_i,\omega_i\}$为第$i$个粒子的空间坐标$\{r_i\}$与自旋坐标$\{\omega_i\}$的集合。根据费米子的泡利原理，交换两粒子的坐标后态矢量变号，即

$$
\begin{equation}
    \ket{x_1\cdots x_i\cdots x_j \cdots x_N} = -\ket{x_1 \cdots x_j \cdots x_i \cdots x_N}
\end{equation}
$$

&emsp;&emsp;根据坐标基矢的定义，只有当两基矢的坐标完全相等时，它们的内积才不为零。又因为对于全同粒子体系，粒子不可分辨，因此对于任意置换$P$，内积都不改变。因此可以将多粒子体系的坐标基矢之间的内积（即坐标基矢在坐标表象下的函数形式）写成如下形式：

$$
\begin{equation}
    \braket{x_1'x_2'\cdots x_N'}{x_1x_2\cdots x_N} = \sum_{P} c_P \prod_{i}\delta(x_i' - x_{P(i)})
\end{equation}
$$

其中求和遍历所有可能的排列。由于我们考虑的是费米子，根据泡利原理，交换任意两个粒子坐标后波函数符号改变，为了构造反对称乘积，可令$c_P = (-1)^P$。因此有

$$
\begin{equation}
    \begin{aligned}
        \braket{x_1'x_2'\cdots x_N'}{x_1x_2\cdots x_N} &= \sum_{P} (-1)^P \prod_{i}\delta(x_i' - x_{P(i)}) \\
        &=
        \begin{vmatrix}
            \delta(x_1'-x_1) & \delta(x_1'-x_2) & \cdots & \delta(x_1'-x_N) \\
            \delta(x_2'-x_1) & \delta(x_2'-x_2) & \cdots & \delta(x_2'-x_N) \\
            \vdots & \vdots & \ddots & \vdots \\
            \delta(x_N'-x_1) & \delta(x_N'-x_2) & \cdots & \delta(x_N'-x_N)
        \end{vmatrix}
    \end{aligned}
\end{equation}
$$

可以得到坐标表象的完备性关系，即

$$
\begin{equation}
    \frac{1}{N!} \int dx_1 \cdots dx_N \ket{x_1x_2\cdots x_N}\bra{x_1x_2\cdots x_N} = \hat{1}
\end{equation}
$$

系数$\frac{1}{N!}$是来源于费米子的不可分辨性。

## 场算符

### 场算符的定义及其对易关系

&emsp;&emsp;引入场算符$\psi^\dagger(x)$，其定义为

$$
\begin{equation}
    \begin{aligned}
        \ket{x_1} &= \psi^\dagger(x_1)\ket{0} \\
        \ket{x_1x_2} &= \psi^\dagger(x_2)\ket{x_1}= \psi^\dagger(x_2)\psi^\dagger(x_1)\ket{0} \\
        \ket{x_1\cdots x_N} &= \psi^\dagger(x_N)\ket{x_1 \cdots x_{N-1}} = \psi^\dagger(x_N)\cdots\psi^\dagger(x_1)\ket{0}
    \end{aligned}
\end{equation}
$$

可以看到，场算符$\psi^\dagger(x)$能将希尔伯特空间$\mathcal{H}_{N}$中的态矢量映射到希尔伯特空间$\mathcal{H}_{N+1}$中，即$\psi^\dagger(x)$的作用是能在$x$处“产生”一个粒子，因此$\psi^\dagger(x)$被称为产生算符（creation operator）。由于我们考虑的是费米子体系，因此体系的态矢量满足交换反对称性：交换体系任意两个粒子，态矢量反号。反映到产生算符上，就要求场算符具有反对易关系，即

$$
\begin{equation}
    \big[\psi^\dagger(x),\psi^\dagger(y)\big]_+ = \psi^\dagger(x)\psi^\dagger(y) + \psi^\dagger(y)\psi^\dagger(x) = 0
\end{equation}
$$

类似地，定义产生算符的共轭为$\psi(x) = \big(\psi^\dagger(x)\big)^\dagger$。将其作用到矢量$\ket{\Psi}$上，与另一矢量$\ket{\Phi}$做内积，可得

$$
\begin{equation}
    \bra{\Phi}\psi(x)\ket{\Psi} = \big(\bra{\Psi}\psi^\dagger(x)\ket{\Phi}\big)^\dagger
\end{equation}
$$

如果$\ket{\Psi} \in \mathcal{H}_{N}$，只有当$\ket{\Phi} \in \mathcal{H}_{N-1}$时，上式右边才可能不为零。因此$\psi(x)$能将希尔伯特空间$\mathcal{H}_{N}$中的态矢量映射到希尔伯特空间$\mathcal{H}_{N-1}$中当$\ket{\Psi}=\ket{0}$时，上式右边的$\psi^\dagger(x)\ket{\Phi} \notin \mathcal{H}_{0}$，因此上式恒等于零，所以$\psi(x)$作用到零矢量后恒等于零，即$\psi(x)$的作用是能在$x$处“湮灭”一个粒子。因此$\psi(x)$被称为湮灭算符（annihilation operator）。同理，由于体系为费米子体系，湮灭算符也满足反对易关系：

$$
\begin{equation}
    \big[\psi(x),\psi(y)\big]_+ = \psi(x)\psi(y) + \psi(y)\psi(x) = 0
\end{equation}
$$

&emsp;&emsp;下面我们来推导产生算符与湮灭算符之间的对易关系。根据上一节的结论（7），根据行列式按行列展开定理，对行列式的第$k$列展开，可得

$$
\begin{equation}
    \begin{aligned}
        \braket{x_1\cdots x_N}{y_1 \cdots y_N}
        &= \bra{x_1\cdots x_{N-1}}\psi(x_N)\ket{y_1\cdots y_N} \\
        &=
        \begin{vmatrix}
            \delta(x_1-y_1) & \delta(x_1-y_2) & \cdots & \delta(x_1-y_N) \\
            \delta(x_2-y_1) & \delta(x_2-y_2) & \cdots & \delta(x_2-y_N) \\
            \vdots & \vdots & \ddots & \vdots \\
            \delta(x_N-y_1) & \delta(x_N-y_2) & \cdots & \delta(x_N-y_N)
        \end{vmatrix} \\
        &= \sum_{k} \delta(x_N - y_k)(-1)^{N+k} \bra{x_1\cdots x_{N-1}}\ket{y_1\cdots y_{k-1} y_{k+1}\cdots y_N} \\
    \end{aligned}
\end{equation}
$$

因此
$$
\begin{equation}
    \begin{aligned}
        \psi(x)\ket{y_1\cdots y_N}
        &= \sum_{k} (-1)^{N+k}\delta(x - y_k)\ket{y_1\cdots y_{k-1} y_{k+1}\cdots y_N}
    \end{aligned}
\end{equation}
$$

$$
\begin{equation}
    \begin{aligned}
        \psi(x)\psi^\dagger(y)\ket{y_1\cdots y_N}
        &= \psi(x) \ket{y_1\cdots y_N y} \\
        &= \delta(x-y)\ket{y_1\cdots y_N} + \sum_{k} (-1)^{N+k}\delta(x - y_k)\ket{y_1\cdots y_{k-1} y_{k+1}\cdots y_N y} \\
        &= \delta(x-y)\ket{y_1\cdots y_N} + \psi^\dagger(y) \sum_{k} (-1)^{N+k}\delta(x - y_k)\ket{y_1\cdots y_{k-1} y_{k+1}\cdots y_N} \\
        &= \delta(x-y)\ket{y_1\cdots y_N} + \psi^\dagger(y) \psi(x) \ket{y_1\cdots y_N}
    \end{aligned}
\end{equation}
$$

由此可得，产生算符与湮灭算符之间的对易关系为

$$
\begin{equation}
    \big[\psi(x),\psi^\dagger(y)\big] = \delta(x-y)
\end{equation}
$$

综上所述，场算符满足如下（反）对易关系：

$$
\begin{equation}
    \boxed{
    \begin{aligned}
        &\big[\psi^\dagger(x),\psi^\dagger(y)\big]_+ = 0 \\
        &\big[\psi(x),\psi(y)\big]_+ = 0 \\
        &\big[\psi(x),\psi^\dagger(y)\big] = \delta(x-y)
    \end{aligned}
}
\end{equation}
$$

### 场算符的表象变换关系

&emsp;&emsp;有了坐标表象下的场算符之后，我们可以通过表象变换将场算符转换到任意表象下。定义离散基$\ket{n}$，将其完备性关系$\hat{1} = \sum\ket{n}\bra{n}$作用到坐标基矢$\ket{x}$上，可得

$$
\begin{equation}
    \begin{aligned}
        \ket{x}
        &= \sum \ket{n}\braket{n}{x}\\
        &= \sum \braket{n}{x} a^\dagger_n \ket{0} \\
        &= \psi^\dagger(x)\ket{0}
    \end{aligned}
\end{equation}
$$

定义$\braket{x}{n}=\varphi_n(x)$为离散基在坐标表象下的波函数，可得

$$
\begin{equation}
    \begin{gathered}
        \psi^\dagger(x) = \sum_{k} \varphi^*_k(x) a^\dagger_k \\
        \psi(x)  = \sum_{k} \varphi_k(x) a_k= (\psi^\dagger(x))^\dagger
    \end{gathered}
\end{equation}
$$

同理，将完备性关系$\hat{1} = \int dx \ket{x}\bra{x}$作用到离散基$\ket{n}$上，可得

$$
\begin{equation}
    \begin{aligned}
        \ket{n}
        &= \int dx \ket{x}\braket{x}{n}\\
        &= \int dx \varphi_n(x) \psi^\dagger(x) \ket{0} \\
        &= a^\dagger_n \ket{0}
    \end{aligned}
\end{equation}
$$

可得

$$
\begin{equation}
    \begin{gathered}
        a^\dagger_n = \int dx \varphi_n(x) \psi^\dagger(x) \\
        a_n = \int dx \varphi^*_n(x) \psi(x)= (a^\dagger_n)^\dagger
    \end{gathered}
\end{equation}
$$

式（19）和式（21）即为连续表象（如坐标表象）的场算符与离散表象（如单电子表象）的场算符之间的转换关系。

通过上述转换关系，我们可以得到离散基基矢在坐标表象下的表达式。利用（21），可得

$$
\begin{equation}
    \begin{aligned}
        \ket{n_1\cdots n_N}
        &= a^\dagger_{n_N} \cdots a^\dagger_{n_1}\ket{0} \\
        &= \int dx_N\cdots dx_1 \varphi_{n_N}(x_N) \psi^\dagger(x_N) \cdots \varphi_{n_1}(x_1) \psi^\dagger(x_1)\ket{0} \\
        &= \int dx_1 \cdots dx_N \varphi_{n_1}(x_1)\cdots \varphi_{n_N}(x_N) \ket{x_1\cdots x_N}
    \end{aligned}
\end{equation}
$$

上式左乘坐标基矢$\bra{x_1\cdots x_N}$，可得

$$
\begin{equation}
    \begin{aligned}
        \braket{x_1\cdots x_N}{n_1\cdots n_N}
        &=  \int dx_1' \cdots dx_N' \varphi_{n_1}(x_1')\cdots \varphi_{n_N}(x_N')  \braket{x_1\cdots x_N}{x_1'\cdots x_N'} \\
        &=  \int dx_1' \cdots dx_N' \varphi_{n_1}(x_1')\cdots \varphi_{n_N}(x_N')  \sum_{P} (-1)^P \prod_{i}\delta(x_i' - x_{P(i)}) \\
        &=
        \begin{vmatrix}
            \varphi_{n_1}(x_1) & \varphi_{n_2}(x_2) & \cdots & \varphi_{n_1}(x_N) \\
            \varphi_{n_2}(x_1) & \varphi_{n_2}(x_2) & \cdots & \varphi_{n_2}(x_N) \\
            \vdots & \vdots & \ddots & \vdots \\
            \varphi_{n_N}(x_1) & \varphi_{n_N}(x_2) & \cdots & \varphi_{n_N}(x_N)
        \end{vmatrix}
    \end{aligned}
\end{equation}
$$

上式即为离散基在坐标表象下的波函数形式。若将$\varphi_n(x)$取为单电子基，则上式即为多电子体系的Slater行列式波函数。

## 用场算符表示哈密顿量

&emsp;&emsp;多电子体系的哈密顿算符在坐标表象下可以表示为

$$
\begin{equation}
    \begin{aligned}
        \hat{H} &= \hat{H}_0 + \hat{H}_{int}\\
        &= \sum_i \hat{h}_i + \frac{1}{2}\sum_{i \neq j} V(x_i,x_j)
    \end{aligned}
\end{equation}
$$

我们可以看到哈密顿量分为两项：单体项和两体项。下面我们来推导这两项的场算符形式。

&emsp;&emsp;定义单电子基$\ket{n}$为单电子哈密顿量$\hat{h}$的本征态，即

$$
\begin{equation}
    \hat{h}\ket{n} = \varepsilon_n \ket{n}
\end{equation}
$$

$$
\begin{equation}
    \begin{aligned}
        \hat{H}_0 \ket{n_1\cdots n_N}
        &= (\varepsilon_{n_1} + \varepsilon_{n_2} + \cdots + \varepsilon_{n_N}) \ket{n_1\cdots n_N}\\
        &= \sum_{i} \varepsilon_{n}\delta_{n n_i} \ket{n_1 \cdots n_N} \\
        &= \sum_{n} \varepsilon_{n} a^\dagger_n a_n \ket{n_1 \cdots n_N}
    \end{aligned}
\end{equation}
$$

因此，单体哈密顿可以表示为

$$
\begin{equation}
    \hat{H}_0 = \sum_{n} \varepsilon_{n} a^\dagger_n a_n
\end{equation}
$$

将式（21）代入上式，可得

$$
\begin{equation}
    \begin{aligned}
        \hat{H}_0
        &= \sum_{n} \varepsilon_{n} \int dx dx' \varphi_n(x) \varphi^*_n(x') \psi^\dagger(x)\psi(x') \\
        &= \int dx dx' \bigg(\sum_{n} \varepsilon_{n} \varphi_n(x) \varphi^*_n(x')\bigg) \psi^\dagger(x)\psi(x') \\
        &= \int dx dx' \bigg(\sum_{n} \braket{x}{n}\varepsilon_{n} \braket{n}{x}\bigg) \psi^\dagger(x)\psi(x') \\
        &= \int dx dx' \bigg(\sum_{n} \braket{x}{n}\hat{h} \braket{n}{x'}\bigg) \psi^\dagger(x)\psi(x') \\
        &= \int dx dx' \bra{x}\hat{h}\ket{x'} \psi^\dagger(x)\psi(x')
    \end{aligned}
\end{equation}
$$

上式的推导过程中用到了单电子基的完备性关系$\hat{1} = \sum \ket{n}\bra{n}$。因此，哈密顿量的单体项可用场算符表示为

$$
\begin{equation}
    \hat{H}_0 = \int dx dx' \bra{x}\hat{h}\ket{x'} \psi^\dagger(x)\psi(x')
\end{equation}
$$

下面来推导两体项。引入密度算符

$$
\begin{equation}
    \hat{n}(x) = \psi^\dagger(x)\psi(x)
\end{equation}
$$

将其作用到坐标基矢上，可得

$$
\begin{equation}
    \begin{aligned}
        \hat{n}(x)\ket{x_1x_2\cdots x_N}
        &=  \psi^\dagger(x)\psi(x)\psi^\dagger(x_N)\cdots\psi^\dagger(x_1)\ket{0}\\
        &= \big(\delta(x-x_N)\psi^\dagger(x) - \psi^\dagger(x_N)\psi^\dagger(x)\psi(x)\big)\psi^\dagger(x_{N-1}) \cdots \psi^\dagger(x_1)\ket{0}\\
        &= \sum_{i} \delta(x-x_i) \ket{x_1x_2\cdots x_N} \\
        &= \rho(x)\ket{x_1x_2\cdots x_N}
    \end{aligned}
\end{equation}
$$

式中$\rho(x) = \sum_i\delta(x-x_i)$为多粒子体系在$x$处的的经典粒子密度。因此算符$\hat{n}(x)$被称为密度算符。

$$
\begin{equation}
    \begin{aligned}
        \hat{H}_{int}\ket{x_1x_2\cdots x_N}
        &= \frac{1}{2}\sum_{i\neq j} V(x_i,x_j)\ket{x_1x_2\cdots x_N} \\
        &= \frac{1}{2}\int dx dx' V(x,x') \sum_{i\neq j} \delta(x-x_i)\delta(x'-x_j)\ket{x_1x_2\cdots x_N} \\
        &= \frac{1}{2}\int dx dx' V(x,x') \big(\hat{n}(x)\hat{n}(x') - \delta(x-x')\hat{n}(x)\big)\ket{x_1x_2\cdots x_N} \\
        &= \frac{1}{2}\int dx dx' V(x,x') \big(\psi^\dagger(x)\psi(x)\psi^\dagger(x')\psi(x') - \delta(x-x')\psi^\dagger(x)\psi(x)\big)\ket{x_1x_2\cdots x_N} \\
        &= \frac{1}{2}\int dx dx' V(x,x') \psi^\dagger(x)\psi^\dagger(x')\psi(x')\psi(x)\ket{x_1x_2\cdots x_N}
    \end{aligned}
\end{equation}
$$

因此

$$
\begin{equation}
    \hat{H}_{int} =  \frac{1}{2}\int dx dx' V(x,x') \psi^\dagger(x)\psi^\dagger(x')\psi(x')\psi(x)
\end{equation}
$$

综上所述，体系的哈密顿量可用场算符表示为

$$
\begin{equation}
    \begin{aligned}
        \hat{H} &= \hat{H}_0 + \hat{H}_{int} \\
        &= \int dxdx' \bra{x}\hat{h}\ket{x'} \psi^\dagger(x)\psi(x') + \frac{1}{2}\int dx dx' V(x,x') \psi^\dagger(x)\psi^\dagger(x')\psi(x')\psi(x)
    \end{aligned}
\end{equation}
$$

对于离散表象，将式（21）代入上式，可得

$$
\begin{equation}
    \hat{H} = \sum_{pq} h_{pq}a^\dagger_p a_q + \frac{1}{2}\sum_{pqrs} \braket{pq|rs} a^\dagger_p a^\dagger_q a_s a_r
\end{equation}
$$

式（34）和式（35）即为哈密顿算符在连续表象和离散表象的场算符表达式

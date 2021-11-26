# Background

The following text is based on the documentation of DDEELS user's guide.

## Generalities and optical response

In the DDA, a particle is described by a set of ``N`` point dipoles interacting with each other. The dipole moment ``\textbf{p}_j`` located at position ``\textbf{r}_j`` is directly proportional to the local
electric field and is given by

```math
\textbf{p}_j = \alpha_j \, \textbf{E}_j^{loc} ,
```

where ``\alpha_j`` is the polarizability of the volume element ``j`` obtained from the bulk dielectric function ``\epsilon_j\left(\omega\right)`` of the material via a Clausius-Mossotti-like relation, and ``\textbf{E}_j^{loc} \equiv \textbf{E}^{loc} \left(\textbf{r}_j\right)`` is the local electric field, i.e. the superposition of the applied field ``\textbf{E}_j^{app}`` and the fields radiating from the ``N-1`` other dipoles. Therefore,

```math
\textbf{p}_j = \alpha_j \left[ \textbf{E}_j^{app} + \sum\limits_{m \neq j} \textbf{A}_{jm} \textbf{p}_{m} \right] ,
```

with the field at the lattice site ``j`` due to a dipole on site ``m`` defined as

```math
\textbf{A}_{jm} \textbf{p}_{m} = -\frac{\exp\left(\mathrm{i}k{r}_{jm}\right)}{{r}_{jm}^{3}} \left\{
    k^2 \textbf{r}_{jm}	\times \left(\textbf{r}_{jm} \times \textbf{p}_{m} \right) +
    \frac{\left(1 - \mathrm{i} k{r}_{jm} \right)}{{r}_{jm}^{2}} \left[
        {r}_{jm}^{2} \textbf{p}_{m} - 3 \textbf{r}_{jm} \left(
            \textbf{r}_{jm}\cdot\textbf{p}_{m}
        \right)
    \right]
\right\},
```

where ``k`` is the wave number and ``{r}_{jm} = | \textbf{r}_j - \textbf{r}_{m} |``. Equation (\ref{eq:dip_j}) is valid for each dipoles ``\textbf{p}_j`` and leads to a system of ``3N`` complex coupled linear equations that can be solved numerically. The electromagnetic response properties can then be computed.

For the investigation of the optical response, an electromagnetic plane wave is considered as applied electric field:

```math
\textbf{E}_j^{app} = \textbf{E}_0 \exp\left( \mathrm{i} \textbf{k} \cdot \textbf{r}_j \right).
```

The extinction cross section is

```math
C_{ext} = \frac{4 \pi k}{|\textbf{E}_{0}|^2} \sum\limits_{j=1}^N \mathrm{Im} \left\{ \textbf{E}_j^{app,*} \cdot \textbf{p}_j \right\}
```

is the optical extinction cross section.


## EELS in the frame of the DDA

To model an EELS experiment, we consider the probe electron as a classical charge ``q`` moving in vacuum with velocity ``\textbf{v}`` parallel to the ``Z``-axis and an impact parameter ``\textbf{r}_{q}`` in the ``XOY``-plan. A dipole at the position ``\textbf{r}_j = \left(x_j, y_j, z_j\right)`` is affected by an applied electric field (in the frequency space) given by

```math
\mathrm{E}_{j,x}^{app} = \frac{q \omega}{2 \pi \epsilon_0 v^2 \gamma}
    \exp\left( \mathrm{i} \omega \frac{z_j}{v} \right)
    \frac{{d}_{j,x}} {{d}_j} {K}_{1}\left(\frac{\omega {d}_j}{\gamma v}\right),
```
```math
\mathrm{E}_{j,y}^{app} = \frac{q \omega}{2 \pi \epsilon_0 v^2 \gamma}
    \exp\left( \mathrm{i} \omega \frac{z_j}{v} \right)
    \frac{{d}_{j,y}} {{d}_j} {K}_{1}\left(\frac{\omega {d}_j}{\gamma v}\right),
```
and
```math
\mathrm{E}_{j,z}^{app} = - \frac{q \omega}{2 \pi \epsilon_0 v^2 \gamma^2}
    \exp\left( \mathrm{i} \omega \frac{z_j}{v} \right)
    \mathrm{i}{K}_{0}\left(\frac{\omega {d}_j}{\gamma v}\right),
```

with ``\gamma = \sqrt{1 - (\frac{v}{c})^2}``, ``\textbf{d}_j = \textbf{r}_j - \textbf{r}_q`` in the ``XOY``-plane, and ``K_m`` is the modified Bessel function of order ``m``. ``\textbf{E}_j^{app}`` depends on the impact parameter (``\textbf{d}_j``).

The total energy ``\left(\Delta E > 0\right)`` exchange between the incident probe charge and the particle (oscillating dipole system) is defined as

```math
\Delta E = \int_0^\infty \hbar \omega \Gamma(\omega) \mathrm{d} (\hbar \omega),
```

where ``\Gamma\left(\omega\right)`` is the loss probability given by

```math
\Gamma_{EELS}(\omega) = \frac{1}{\pi \hbar^2} \sum\limits_{j=1}^{N} \mathrm{Im} \left\{ \textbf{E}_j^{app,*} \cdot \textbf{p}_j \right\}.
```

In the case of one isolated dipole (``N=1``), we have ``\textbf{p}_{1} = \alpha_1 \textbf{E}_{1}^{app}`` and

```math
\Gamma_{EELS}(\omega) = \frac{1}{\pi \hbar^2} \mathrm{Im} \left\{ \alpha_1 \right\} \left| \textbf{E}_{1}^{app} \right|^2,
```

as found by Echenique et al. for the dipolar loss of a dielectric sphere.


## Cathodoluminescence in the frame of the DDA

For cathodoluminescence, light emission probability is given by (Eqs. (\ref{Ex}-\ref{lossener}) remain valid)

```math
\Gamma_{CL}(\omega) = \frac{c \, r_s^2}{4 \pi^2 \hbar^2 \omega} \int \mathrm{d} \Omega | \mathbf{E}_{ind}(\mathbf{r}_s,\omega)|
```

where the integral is performed on the far-field induced field, i.e., on a sphere with ``kr_s >> 1``. The induced field is obtained by adding the field from each dipole moments (``\mathbf{r}_{sj}= \mathbf{r}_s - \mathbf{r}_j`` and Eq. (\ref{ajm}))

```math
\mathbf{E}_{ind}(\mathbf{r}_s,\omega)=\sum_j \mathbf{A}_{sj} \mathbf{p}_j
```


## Substrate effect

The influence of a semi-infinite substrate is taken into account using the classical image model. In our EELS description, the expression of the probe charge electric field (Eqs. (\ref{Ex}-\ref{Ez})) requires an electron trajectory in vacuum.

First, let the substrate surface be the ``YOZ``-plan (``x=0``). The polarization field of the semi-infinite substrate due to a charged particle at ``(x_q(t),y_q(t),z_q(t))`` can be described by an image-charge ``{q}^{'}`` located in the substrate at ``(-x_q(t),y_q(t),z_q(t))``. The total applied field at the lattice site ``j`` becomes the sum of the field of the actual charge and its image:

```math
\textbf{E}_j^{app} = \textbf{E}_j^{app}[q] + \textbf{E}_j^{app}[q'] ,
```

where

```math
{q}' = \frac{1-\epsilon_{sub}}{1+\epsilon_{sub}} \, q ,
```

with ``\epsilon_{sub}`` the dielectric function of the substrate.

Second, the substrate response to the dipole moments is also considered via the image method. Equation (\ref{eq:dip_j}) becomes

```math
\textbf{p}_j = \alpha_j \left[ \textbf{E}_j^{app} + \sum\limits_{m \neq j} \textbf{A}_{jm} \textbf{p}_{m} + \sum_{m} \textbf{A}_{jm}' \textbf{p}_{m}' \right] ,
```

where the last term represents the image-dipole contribution to the local field and with

```math
\mathrm{p}_{j,x}^{'} = -\frac{1 - \epsilon_{sub}}{1 + \epsilon_{sub}} \, \mathrm{p}_{j,x},
```
```math
\mathrm{p}_{j,y}^{'} =  \frac{1 - \epsilon_{sub}}{1 + \epsilon_{sub}} \, \mathrm{p}_{j,y}
```
and
```math
\mathrm{p}_{j,z}^{'} =  \frac{1 - \epsilon_{sub}}{1 + \epsilon_{sub}} \, \mathrm{p}_{j,z},
```

where the image dipole is located at ``\textbf{r}_j^{'} = \left(-x_j, y_j, z_j\right)`` for a dipole moment ``\textbf{p}_j`` at ``{\textbf{r}_j} = \left(x_j, y_j, z_j\right)``.

Once the linear system is solved, both ``\textbf{p}_j`` (nanoparticle polarization) and ``\mathrm{p}_j^{'}`` (substrate polarization) are known, and both contributions must be included in the loss probability (Eq. (\ref{eq:loss})).

We emphasize here that the present formulation (and, in particular, Eq. (\ref{chfield})) is only valid when the electron does not penetrate into the substrate.

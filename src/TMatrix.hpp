/**
 * @file TMatrix.hpp
 *
 * @brief T-matrix class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TMATRIX_HPP
#define TMATRIX_HPP

#include "Configuration.hpp"
#include "Matrix.hpp"
#include "SpecialFunctions.hpp"

#include <complex>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief T-matrix and auxiliary variables required to compute it.
 *
 * The T-matrix formalism is based on rewriting the incoming @f$(i)@f$, local
 * @f$(l)@f$ and scattered @f$(s)@f$ electromagnetic waves (photons) in terms of
 * spherical wave function expansions, with spherical basis functions
 * @f$Rg\vec{M}_{mn}, Rg\vec{N}_{mn}, \vec{M}_{mn}@f$ and @f$\vec{N}_{mn}@f$
 * (introduced below; the @f$Rg@f$ basis functions are the real variants of the
 * general basis functions):
 * @f[
 *    \vec{E}^{(i)}(r, \theta{}, \phi{}) = \sum_{n=1}^{n_{max}} \sum_{m=-n}^{n}
 *        \left( a_{mn} Rg\vec{M}_{mn}(kr, \theta{}, \phi{}) +
 *               b_{mn} Rg\vec{N}_{mn}(kr, \theta{}, \phi{}) \right),
 * @f]
 * @f[
 *    \vec{E}^{(l)}(r, \theta{}, \phi{}) = \sum_{n=1}^{n_{max}} \sum_{m=-n}^{n}
 *        \left( c_{mn} Rg\vec{M}_{mn}(km_rr, \theta{}, \phi{}) +
 *               d_{mn} Rg\vec{N}_{mn}(km_rr, \theta{}, \phi{}) \right),
 * @f]
 * and
 * @f[
 *     \vec{E}^{(s)}(r, \theta{}, \phi{}) = \sum_{n=1}^{n_{max}} \sum_{m=-n}^{n}
 *        \left( p_{mn} \vec{M}_{mn}(kr, \theta{}, \phi{}) +
 *               q_{mn} \vec{N}_{mn}(kr, \theta{}, \phi{}) \right).
 * @f]
 * In these expressions, all spatial coordinates are expressed in spherical
 * coordinates, @f$r, \theta{}, \phi{}@f$, and the spherical wave function
 * expansions are assumed to converge for some finite order
 * @f$n_{max} \geq{} 1@f$. @f$k = \frac{2\pi{}}{\lambda{}}@f$ is the wave number
 * of the incident and outgoing radiation (assuming Rayleigh scattering), and
 * @f$m_r@f$ is the refractive index of the dust grain that causes the
 * scattering.
 *
 * The spherical basis functions for the expansion differ for different authors.
 * We follow Mishchenko, Travis & Mackowski, 1996, J. Quant. Spectrosc.
 * Radiat. Transfer, 55, 535 (https://doi.org/10.1016/0022-4073(96)00002-7) and
 * use the following expressions:
 * @f[
 *    \vec{M}_{mn}(kr, \theta{}, \phi{}) =
 *      (-1)^m \frac{1}{\sqrt{4\pi{}}} \sqrt{\frac{2n+1}{n(n+1)}} e^{im\phi{}}
 *      h^{(1)}_n(kr) \left[
 *        \frac{im}{\sin(\theta)} d^n_{0m}(\theta{}) \hat{\theta{}}
 *        - \frac{d}{d\theta{}} d^n_{0m}(\theta{}) \hat{\phi{}}
 *      \right] \\ = \vec{\nabla{}} \times{} \left(
 *        (-1)^m \frac{1}{\sqrt{4\pi{}}} \sqrt{\frac{2n+1}{n(n+1)}}
 *        r h^{(1)}_n(kr) e^{im\phi{}} d^n_{0m}(\theta{}) \hat{r}
 *      \right),
 * @f]
 * and
 * @f[
 *    \vec{N}_{mn}(kr, \theta{}, \phi{}) =
 *      (-1)^m \frac{1}{\sqrt{4\pi{}}} \sqrt{\frac{2n+1}{n(n+1)}} e^{im\phi{}}
 *      \left[
 *        \frac{n(n+1)}{kr} h^{(1)}_n(kr) d^n_{0m}(\theta{}) \hat{r}
 *        + \frac{[kr h^{(1)}_n(kr)]'}{kr} \left(
 *          \frac{d}{d\theta{}} d^n_{0m}(\theta{}) \hat{\theta{}}
 *          + \frac{im}{\sin(\theta{})} d^n_{0m}(\theta{}) \hat{\phi{}}
 *        \right)
 *      \right] \\ = \frac{1}{k} \vec{\nabla{}} \times{} \left(
 *        \vec{\nabla{}} \times{} \left(
 *          (-1)^m \frac{1}{\sqrt{4\pi{}}} \sqrt{\frac{2n+1}{n(n+1)}}
 *          r h^{(1)}_n(kr) e^{im\phi{}} d^n_{0m}(\theta{}) \hat{r}
 *        \right)
 *      \right).
 * @f]
 * In these expressions, @f$\hat{r}, \hat{\theta{}}@f$ and @f$\hat{\phi{}}@f$
 * represent the unit vectors in spherical coordinates, and the curl operator
 * is given by
 * (https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates#Del_formula)
 * @f[
 *    \vec{\nabla{}} \times{} \left(
 *      A_r \hat{r} + A_\theta{} \hat{\theta{}} + A_\phi{} \hat{\phi{}}
 *    \right) =
 *      \frac{1}{r\sin(\theta{})} \left(
 *        \frac{\partial{}}{\partial{}\theta{}} \left(
 *          A_\phi{} \sin(\theta{})
 *        \right) - \frac{\partial{}A_\theta{}}{\partial{}\phi{}}
 *      \right) \hat{r} + \frac{1}{r} \left(
 *        \frac{1}{\sin(\theta{})} \frac{\partial{}A_r}{\partial{}\phi{}}
 *        - \frac{\partial{}}{\partial{}r} \left( r A_\phi{} \right)
 *      \right) \hat{\theta{}} + \frac{1}{r} \left(
 *        \frac{\partial{}}{\partial{}r}\left( r A_\theta{} \right)
 *        - \frac{\partial{}A_r}{\partial{}\theta{}}
 *      \right) \hat{\phi{}}.
 * @f]
 *
 * @f$h^{(1)}_n(x)@f$ are the spherical Hankel functions of the first kind,
 * defined as
 * @f[
 *    h^{(1)}_n(x) = j_n(x) + i y_n(x),
 * @f]
 * where @f$j_n(x)@f$ and @f$y_n(x)@f$ are the spherical Bessel functions of
 * respectively the first and second kind, and @f$d^n_{0m}(\theta{})@f$ are
 * the Wigner @f$d@f$-functions (see SpecialFunctions::wigner_dn_0m()) that
 * satisfy
 * @f[
 *    \sin^2(\theta{}) \frac{d^2}{d\cos(\theta{})^2} d^n_{0m}(\theta{})
 *    - 2\cos(\theta{}) \frac{d}{d\cos(\theta{})} d^n_{0m}(\theta{})
 *    + \left(
 *      n(n+1) - \frac{m^2}{\sin^2(\theta{})}
 *    \right) d^n_{0m}(\theta{}) = 0,
 * @f]
 * a variant of the general Legendre equation
 * (https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).
 *
 * Note that the spherical basis functions @f$\vec{M}_{mn}@f$ and
 * @f$\vec{N}_{mn}@f$ just like the electromagnetic fields satisfy the
 * vector Helmholtz equation
 * @f[
 *    \nabla{}^2\vec{A} + k^2\vec{A} = 0,
 * @f]
 * which also leads to
 * (https://en.wikipedia.org/wiki/Vector_calculus_identities#Second_derivatives)
 * @f[
 *    \vec{\nabla{}} \times{} \left(
 *      \vec{\nabla{}} \times{} \vec{M}_{mn}
 *    \right) = \vec{\nabla{}}.\left(\vec{\nabla{}}\times{}\vec{M}_{mn}\right)
 *              -\nabla{}^2 \vec{M}_{mn} = k^2 \vec{M}_{mn},
 * @f]
 * from which we can see that
 * @f[
 *    \frac{1}{k} \vec{\nabla{}} \times{} \vec{N}_{mn} = \vec{M}_{mn},
 * @f]
 * while we already had
 * @f[
 *    \frac{1}{k} \vec{\nabla{}} \times{} \vec{M}_{mn} = \vec{N}_{mn}.
 * @f]
 *
 * The expressions for @f$Rg\vec{M}_{mn}@f$ and @f$Rg\vec{N}_{mn}@f$ are the
 * same as for @f$\vec{M}_{mn}@f$ and @f$\vec{N}_{mn}@f$, except that all
 * occurences of the spherical Hankel functions of the first kind are replaced
 * by spherical Bessel functions of the first kind (as the incoming field is
 * assumed to be real). Note that the coordinate @f$m_rkr@f$ that appears in
 * the expansion for the internal field can be a complex number (since @f$m_r@f$
 * is generally a complex number), so that we also require spherical Bessel
 * functions of the first kind for complex arguments.
 *
 * The expansion coefficients @f$a_{mn}@f$ and @f$b_{mn}@f$ are assumed to be
 * known, while the expansion coefficients @f$c_{mn}, d_{mn}, p_{mn}@f$ and
 * @f$q_{mn}@f$ are not. The T-matrix, or transition matrix, is a matrix that
 * links the expansion coefficients of the incoming field to those of the
 * scattered field:
 * @f[
 *    p_{mn} = \sum_{n'=1}^{n_{max}} \sum_{m'=-n'}^{n'}
 *      \left(
 *        T^{(11)}_{mnm'n'} a_{m'n'} + T^{(12)}_{mnm'n'} b_{m'n'}
 *      \right),
 * @f]
 * and
 * @f[
 *    q_{mn} = \sum_{n'=1}^{n_{max}} \sum_{m'=-n'}^{n'}
 *      \left(
 *        T^{(21)}_{mnm'n'} a_{m'n'} + T^{(22)}_{mnm'n'} b_{m'n'}
 *      \right).
 * @f]
 * The T-matrix only depends on the wavelength of the radiation and the
 * material properties (refractive index) and geometry of the dust. It can be
 * computed numerically by assuming that the expansions for the incoming and
 * scattered waves are valid far from the dust particle, and by assuming that
 * the Maxwell equations hold throughout the scattering process. This technique
 * is called the Extended Boundary Condition Method (EBCM), and was introduced
 * by Waterman, 1971, Physical Review D, 3, 825
 * (https://doi.org/10.1103/PhysRevD.3.825).
 *
 * Skipping the detailed derivation, it can be shown that the transition matrix
 * can be computed from
 * @f[
 *    T = -RgQ Q^{-1},
 * @f]
 * where the matrix Q is given by (ignoring constant prefactors that do not
 * contribute to @f$T@f$)
 * @f[
 *    Q_{mnm'n'} = k \int{} d\vec{\sigma{}} . \left[
 *        \left(
 *          \vec{\nabla{}} \times{} Rg\vec{X}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \right) \times{} \vec{X}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{X}_{m'n'}(km_rr, \theta{}, \phi{}) \times{} \left(
 *          \vec{\nabla{}} \times{} \vec{X}_{mn}(kr, \theta{}, \phi{})
 *        \right)
 *    \right],
 * @f]
 * with @f$\vec{X}_{mn} = \vec{M}_{mn}@f$ or @f$\vec{N}_{mn}@f$ depending on the
 * specific component of the matrix under consideration, and the @f$RgQ@f$
 * matrix is given by the same expression where @f$M_{mn}@f$ and @f$N_{mn}@f$
 * are replaced with @f$RgM_{mn}@f$ and @f$RgN_{mn}@f$ respectively. The surface
 * element @f$d\vec{\sigma{}}@f$ is given by
 * @f[
 *    d\vec{\sigma{}} = r^2 \sin(\theta{}) d\theta{} d\phi{} \left(
 *      \hat{r}
 *      - \frac{1}{r}\frac{\partial{}}{\partial{}\theta{}}
 *          r(\theta{}, \phi{}) \hat{\theta{}}
 *      - \sin(\theta{}) \frac{1}{r}\frac{\partial{}}{\partial{}\phi{}}
 *          r(\theta{}, \phi{}) \hat{\phi{}}
 *    \right),
 * @f]
 * which for a spheroid (@f$r(\theta{},\phi{}) = r(\theta{})@f$) reduces to
 * @f[
 *    d\vec{\sigma{}} = r^2 \sin(\theta{}) d\theta{} d\phi{} \left(
 *      \hat{r} - \frac{1}{r}\frac{d}{d\theta{}} r(\theta{}) \hat{\theta{}}
 *    \right).
 * @f]
 *
 * To compute the matrix @f$Q@f$, it makes sense to decompose it into four
 * quarters:
 * @f[
 *    Q = \begin{pmatrix}
 *      P & R \\
 *      S & U
 *    \end{pmatrix},
 * @f]
 * with
 * @f[
 *    P = k \int{} d\vec{\sigma{}} . \left[
 *        \left(
 *          \vec{\nabla{}} \times{} Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \right) \times{} \vec{M}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{}) \times{} \left(
 *          \vec{\nabla{}} \times{} \vec{M}_{mn}(kr, \theta{}, \phi{})
 *        \right)
 *    \right] \\ = k^2 \int{} d\vec{\sigma{}} . \left[
 *      Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \times{} \vec{M}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{}) \times{}
 *          \vec{N}_{mn}(kr, \theta{}, \phi{})
 *    \right],
 * @f]
 * @f[
 *    R = k \int{} d\vec{\sigma{}} . \left[
 *        \left(
 *          \vec{\nabla{}} \times{} Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \right) \times{} \vec{M}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{}) \times{} \left(
 *          \vec{\nabla{}} \times{} \vec{M}_{mn}(kr, \theta{}, \phi{})
 *        \right)
 *      \right] \\ = k^2 \int{} d\vec{\sigma{}} . \left[
 *        Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \times{} \vec{M}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{}) \times{}
 *          \vec{N}_{mn}(kr, \theta{}, \phi{})
 *      \right],
 * @f]
 * @f[
 *    S = k \int{} d\vec{\sigma{}} . \left[
 *        \left(
 *          \vec{\nabla{}} \times{} Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \right) \times{} \vec{N}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{}) \times{} \left(
 *          \vec{\nabla{}} \times{} \vec{N}_{mn}(kr, \theta{}, \phi{})
 *        \right)
 *      \right] \\ = k^2 \int{} d\vec{\sigma{}} . \left[
 *        Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \times{} \vec{N}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{}) \times{}
 *        \vec{M}_{mn}(kr, \theta{}, \phi{})
 *    \right],
 * @f]
 * and
 * @f[
 *    U = k \int{} d\vec{\sigma{}} . \left[
 *        \left(
 *          \vec{\nabla{}} \times{} Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \right) \times{} \vec{N}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{}) \times{} \left(
 *          \vec{\nabla{}} \times{} \vec{N}_{mn}(kr, \theta{}, \phi{})
 *        \right)
 *      \right] \\ = k^2 \int{} d\vec{\sigma{}} . \left[
 *        Rg\vec{M}_{m'n'}(km_rr, \theta{}, \phi{})
 *        \times{} \vec{N}_{mn}(kr, \theta{}, \phi{})
 *        - Rg\vec{N}_{m'n'}(km_rr, \theta{}, \phi{}) \times{}
 *          \vec{M}_{mn}(kr, \theta{}, \phi{})
 *      \right].
 * @f]
 *
 * We have
 * @f[
 *    d\vec{\sigma{}}.\left(
 *      \vec{M}_{mn}(k_1r) \times{} \vec{M}_{m'n'}(k_2r)
 *    \right) =
 *      (-1)^{m+m'} \frac{1}{4\pi{}}
 *      \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
 *      e^{i(m+m')\phi{}} h_n^{(1)}(k_1r) h_{n'}^{(1)}(k_2r)
 *      \left(d\vec{\sigma{}}.\hat{r}\right) \left(
 *        -\frac{im}{\sin(\theta{})} d^n_{0m}(\theta{})
 *          \frac{d}{d\theta{}} d^{n'}_{0m'}(\theta{})
 *        +\frac{im'}{\sin(\theta{})} d^{n'}_{0m'}(\theta{})
 *          \frac{d}{d\theta{}} d^n_{0m}(\theta{})
 *      \right) \\ =
 *      (-1)^{m} \frac{1}{4\pi{}}
 *      \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
 *      e^{i(m-m')\phi{}} h_n^{(1)}(k_1r) h_{n'}^{(1)}(k_2r)
 *      \left(d\vec{\sigma{}}.\hat{r}\right) \left(
 *        -\frac{im}{\sin(\theta{})} d^n_{0m}(\theta{})
 *          \frac{d}{d\theta{}} d^{n'}_{0m'}(\theta{})
 *        -\frac{im'}{\sin(\theta{})} d^{n'}_{0m'}(\theta{})
 *          \frac{d}{d\theta{}} d^n_{0m}(\theta{})
 *      \right),
 * @f]
 * where for the last step we used
 * @f[
 *    d^n_{0m}(\theta{}) = (-1)^{-m} d^n_{0-m}(\theta{}).
 * @f]
 * Similarly, we have
 * @f[
 *    d\vec{\sigma{}}.\left(
 *      \vec{M}_{mn}(k_1r) \times{} \vec{N}_{m'n'}(k_2r)
 *    \right) \\=
 *      (-1)^{m+m'} \frac{1}{4\pi{}}
 *      \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
 *      e^{i(m+m')\phi{}} d\vec{\sigma{}}. \left[
 *        h_n^{(1)}(k_1r) \frac{[k_2rh_{n'}^{(1)}(k_2r)]'}{k_2r} \left(
 *         \frac{-mm'}{\sin^2(\theta{})} d^n_{0m}(\theta{})
 *           d^{n'}_{0m'}(\theta{})
 *        + \frac{d}{d\theta{}}d^n_{0m}(\theta{})
 *          \frac{d}{d\theta{}} d^{n'}_{0m'}(\theta{})
 *      \right) \hat{r} - h_n^{(1)}(k_1r)
 *        \frac{d}{d\theta{}}d^n_{0m}(\theta{})
 *        \frac{n'(n'+1)}{k_2r} h_{n'}^{(1)}(k_2r) d^{n'}_{0m'}(\theta{})
 *        \hat{\theta{}}
 *      \right] \\=
 *      (-1)^{m} \frac{1}{4\pi{}}
 *      \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
 *      e^{i(m-m')\phi{}} d\vec{\sigma{}}. \left[
 *        h_n^{(1)}(k_1r) \frac{[k_2rh_{n'}^{(1)}(k_2r)]'}{k_2r} \left(
 *         \frac{mm'}{\sin^2(\theta{})} d^n_{0m}(\theta{})
 *           d^{n'}_{0m'}(\theta{})
 *        + \frac{d}{d\theta{}}d^n_{0m}(\theta{})
 *          \frac{d}{d\theta{}} d^{n'}_{0m'}(\theta{})
 *      \right) \hat{r} - h_n^{(1)}(k_1r)
 *        \frac{d}{d\theta{}}d^n_{0m}(\theta{})
 *        \frac{n'(n'+1)}{k_2r} h_{n'}^{(1)}(k_2r) d^{n'}_{0m'}(\theta{})
 *        \hat{\theta{}}
 *      \right],
 * @f]
 * where we have ignored the @f$\hat{\phi{}}@f$ factor that drops out in the
 * inner product with the surface element, and
 * @f[
 *    d\vec{\sigma{}}.\left(
 *      \vec{N}_{mn}(k_1r) \times{} \vec{N}_{m'n'}(k_2r)
 *    \right) = (-1)^{m+m'} \frac{1}{4\pi{}}
 *      \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
 *      e^{i(m+m')\phi{}} d\vec{\sigma{}}. \left[
 *        \frac{[k_1rh_{n}^{(1)}(k_1r)]'}{k_1r}
 *        \frac{[k_2rh_{n'}^{(1)}(k_2r)]'}{k_2r} \left(
 *          \frac{d}{d\theta{}} d^n_{0m}(\theta{})
 *            \frac{im'}{\sin(\theta{})} d^{n'}_{0m'}(\theta{})
 *          - \frac{im}{\sin(\theta{})} d^{n}_{0m}(\theta{})
 *            \frac{d}{d\theta{}} d^{n'}_{0m'}(\theta{})
 *        \right) \hat{r} \\+ d^{n}_{0m}(\theta{}) d^{n'}_{0m'}(\theta{}) \left(
 *          \frac{[k_1rh_{n}^{(1)}(k_1r)]'}{k_1r} \frac{im}{\sin(\theta{})}
 *            \frac{n'(n'+1)}{k_2r} h_{n'}^{(1)}(k_2r)
 *          - \frac{n(n+1)}{k_1r} h_{n}^{(1)}(k_1r)
 *            \frac{[k_2rh_{n'}^{(1)}(k_2r)]'}{k_2r} \frac{im'}{\sin(\theta{})}
 *        \right) \hat{\theta{}}
 *      \right] \\ = (-1)^{m} \frac{1}{4\pi{}}
 *      \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
 *      e^{i(m-m')\phi{}} d\vec{\sigma{}}. \left[
 *        \frac{[k_1rh_{n}^{(1)}(k_1r)]'}{k_1r}
 *        \frac{[k_2rh_{n'}^{(1)}(k_2r)]'}{k_2r} \left(
 *          - \frac{d}{d\theta{}} d^n_{0m}(\theta{})
 *            \frac{im'}{\sin(\theta{})} d^{n'}_{0m'}(\theta{})
 *          - \frac{im}{\sin(\theta{})} d^{n}_{0m}(\theta{})
 *            \frac{d}{d\theta{}} d^{n'}_{0m'}(\theta{})
 *        \right) \hat{r} \\+ d^{n}_{0m}(\theta{}) d^{n'}_{0m'}(\theta{}) \left(
 *          \frac{[k_1rh_{n}^{(1)}(k_1r)]'}{k_1r} \frac{im}{\sin(\theta{})}
 *            \frac{n'(n'+1)}{k_2r} h_{n'}^{(1)}(k_2r)
 *          + \frac{n(n+1)}{k_1r} h_{n}^{(1)}(k_1r)
 *            \frac{[k_2rh_{n'}^{(1)}(k_2r)]'}{k_2r} \frac{im'}{\sin(\theta{})}
 *        \right) \hat{\theta{}}
 *      \right]
 * @f]
 *
 * When integrating out these expressions over all angles @f$\phi{}@f$, we
 * can make use of the identity
 * @f[
 *    \delta{}_{m,m'} =
 *      \frac{1}{2\pi{}} \int_0^{2\pi{}} e^{i(m-m')\phi{}} d\phi{}
 * @f]
 * to get rid of all components for which @f$m\neq{}m'@f$.
 */
class TMatrix {
private:
  /*! @brief Maximum order of the spherical basis functions, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of components in one quarter of the T-matrix,
   *  @f$L_{max} = n_{max} (n_{max} + 2) @f$. */
  const uint_fast32_t _Lmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Precomputed factors @f$n(n+1)@f$ (array of size @f$n_{max}@f$). */
  std::vector<float_type> _an;

  /*! @brief Precomputed factors @f$\sqrt{\frac{2n+1}{n(n+1)}}@f$ (array of size
   *  @f$n_{max}@f$). */
  std::vector<float_type> _dd;

  /*! @brief Precomputed factors
   *  @f$\frac{1}{2}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$
   *  (@f$n_{max}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _ann;

  /*! @brief Precomputed factors @f$\cos(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _costheta;

  /*! @brief Precomputed factors @f$\frac{1}{\sin(\theta{})}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _sinthetainv;

  /*! @brief Precomputed factors @f$\frac{1}{\sin^2(\theta{})}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _sintheta2inv;

  /*! @brief Gauss-Legendre weights for the roots @f$\cos(\theta{})@f$ (array of
   *  size @f$2n_{GL}@f$). */
  std::vector<float_type> _weights;

  /*! @brief Precomputed factors @f$r^2(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _r2;

  /*! @brief Precomputed factors @f$\frac{1}{r(\theta{})}\frac{d}{d\theta{}}
   *  r(\theta{})@f$  (array of size @f$2n_{GL}@f$). */
  std::vector<float_type> _dr_over_r;

  /*! @brief Precomputed factors @f$kr@f$ (array of size @f$2n_{GL}@f$). */
  std::vector<float_type> _kr;

  /*! @brief Precomputed factors @f$\frac{1}{kr}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _krinv;

  /*! @brief Precomputed factors @f$km_rr@f$ (array of size @f$2n_{GL}@f$). */
  std::vector<std::complex<float_type>> _krmr;

  /*! @brief Precomputed factors @f$\frac{1}{km_rr}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<std::complex<float_type>> _krmrinv;

  /*! @brief Wavenumber, @f$k = \frac{2\pi{}}{\lambda{}}@f$. */
  const float_type _k;

  /*! @brief Wavenumber squared, @f$k^2@f$. */
  const float_type _k2;

  /*! @brief Wavenumber squared times refractive index, @f$m_rk^2@f$. */
  const std::complex<float_type> _k2mr;

  /*! @brief Bessel functions @f$j_n(kr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<float_type> _jkr;

  /*! @brief Bessel functions @f$y_n(kr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<float_type> _ykr;

  /*! @brief Bessel function derivatives @f$\frac{[krj_n(kr)]'}{kr}@f$
   *  (@f$2n_{GL}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _djkr;

  /*! @brief Bessel function derivatives @f$\frac{[kry_n(kr)]'}{kr}@f$
   *  (@f$2n_{GL}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _dykr;

  /*! @brief Bessel functions @f$j_n(km_rr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<std::complex<float_type>> _jkrmr;

  /*! @brief Bessel function derivatives
   *  @f$\frac{[km_rrj(km_rr)]'}{km_rr}@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<std::complex<float_type>> _djkrmr;

  /*! @brief T-matrix itself
   *  (@f$2n_{max}(n_{max}+2)\times{}2n_{max}(n_{max}+2)@f$ matrix). */
  Matrix<std::complex<float_type>> _T;

public:
  /**
   * @brief T-matrix.
   *
   * @param wavelength Wavelength of incident radiation, @f$\lambda{}@f$.
   * @param refractive_index Refractive index of the material, @f$m_r@f$.
   * @param R_V Equal volume sphere radius, @f$R_V@f$.
   * @param axis_ratio Axis ratio of the spheroid, @f$d = \frac{a}{b}@f$.
   * @param nmax Maximum order of spherical basis, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline TMatrix(const float_type wavelength,
                 const std::complex<float_type> refractive_index,
                 const float_type R_V, const float_type axis_ratio,
                 const uint_fast32_t nmax, const uint_fast32_t ngauss)
      : _nmax(nmax), _Lmax(nmax * (nmax + 2)), _ngauss(ngauss),
        _an(nmax, float_type(0.)), _dd(nmax, float_type(0.)), _ann(nmax, nmax),
        _costheta(2 * ngauss, float_type(0.)),
        _sinthetainv(2 * ngauss, float_type(0.)),
        _sintheta2inv(2 * ngauss, float_type(0.)),
        _weights(2 * ngauss, float_type(0.)), _r2(2 * ngauss, float_type(0.)),
        _dr_over_r(2 * ngauss, float_type(0.)), _kr(2 * ngauss, float_type(0.)),
        _krinv(2 * ngauss, float_type(0.)), _krmr(2 * ngauss, float_type(0.)),
        _krmrinv(2 * ngauss, float_type(0.)), _k(2. * M_PI / wavelength),
        _k2(_k * _k), _k2mr(refractive_index * _k2), _jkr(2 * ngauss, nmax),
        _ykr(2 * ngauss, nmax), _djkr(2 * ngauss, nmax),
        _dykr(2 * ngauss, nmax), _jkrmr(2 * ngauss, nmax),
        _djkrmr(2 * ngauss, nmax), _T(2 * _Lmax, 2 * _Lmax) {

    for (uint_fast32_t ni = 0; ni < nmax; ++ni) {
      const float_type nn((ni + 2.) * (ni + 1.));
      _an[ni] = nn;
      const float_type d = sqrt((2. * (ni + 1.) + 1.) / nn);
      _dd[ni] = d;
      for (uint_fast32_t nj = 0; nj < ni + 1; ++nj) {
        const float_type ddd = 0.5 * d * _dd[nj];
        _ann(ni, nj) = ddd;
        _ann(nj, ni) = ddd;
      }
    }
    SpecialFunctions::get_gauss_legendre_points_and_weights(
        2 * ngauss, _costheta, _weights);
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      const float_type this_sintheta2inv =
          1. / (1. - _costheta[ig] * _costheta[ig]);
      _sintheta2inv[ig] = this_sintheta2inv;
      _sintheta2inv[2 * ngauss - ig - 1] = this_sintheta2inv;
      const float_type this_sinthetainv = sqrt(this_sintheta2inv);
      _sinthetainv[ig] = this_sinthetainv;
      _sinthetainv[2 * ngauss - ig - 1] = this_sinthetainv;
    }
    SpecialFunctions::get_r_dr_spheroid(_costheta, R_V, axis_ratio, _r2,
                                        _dr_over_r);
    const std::complex<float_type> mrinv = float_type(1.) / refractive_index;
    for (uint_fast32_t i = 0; i < 2 * ngauss; ++i) {
      const float_type r = sqrt(_r2[i]);
      _kr[i] = _k * r;
      _krmr[i] = refractive_index * _kr[i];
      _krinv[i] = 1. / _kr[i];
      _krmrinv[i] = mrinv * _krinv[i];
    }
    for (uint_fast32_t ig = 0; ig < 2 * ngauss; ++ig) {
      SpecialFunctions::spherical_j_jdj_array(nmax, _kr[ig], _jkr.get_row(ig),
                                              _djkr.get_row(ig));
      SpecialFunctions::spherical_y_ydy_array(nmax, _kr[ig], _ykr.get_row(ig),
                                              _dykr.get_row(ig));
      SpecialFunctions::spherical_j_jdj_array(
          nmax, _krmr[ig], _jkrmr.get_row(ig), _djkrmr.get_row(ig));
    }

    const uint_fast32_t nmax2 = 2 * nmax;

    std::vector<int_fast8_t> signs(nmax2);
    Matrix<float_type> wigner_d(2 * ngauss, nmax);
    Matrix<float_type> dwigner_d(2 * ngauss, nmax);
    std::vector<float_type> wr2(ngauss);
    Matrix<std::complex<float_type>> J12(nmax, nmax);
    Matrix<std::complex<float_type>> J21(nmax, nmax);
    Matrix<std::complex<float_type>> RgJ12(nmax, nmax);
    Matrix<std::complex<float_type>> RgJ21(nmax, nmax);
    Matrix<std::complex<float_type>> Q(nmax2, nmax2);
    Matrix<std::complex<float_type>> RgQ(nmax2, nmax2);

    int_fast8_t si = 1;
    for (uint_fast32_t m = 0; m < nmax2; ++m) {
      si = -si;
      signs[m] = si;
    }
    for (uint_fast32_t ig = 1; ig < ngauss + 1; ++ig) {
      const uint_fast32_t i1 = ngauss + ig;
      const uint_fast32_t i2 = ngauss - ig + 1;
      std::vector<float_type> dv1(nmax), dv2(nmax);
      SpecialFunctions::wigner_dn_0m(_costheta[i1 - 1], nmax, 0, &dv1[0],
                                     &dv2[0]);
      for (uint_fast32_t n = 0; n < nmax; ++n) {
        si = signs[n];
        wigner_d(i1 - 1, n) = dv1[n];
        wigner_d(i2 - 1, n) = si * dv1[n];
        dwigner_d(i1 - 1, n) = dv2[n];
        dwigner_d(i2 - 1, n) = -si * dv2[n];
      }
    }
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      wr2[ig] = _weights[ig] * _r2[ig];
    }
    for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
      // n1 * (n1 + 1)
      const float_type n1n1p1 = _an[n1 - 1];
      for (uint_fast32_t n2 = 1; n2 < nmax + 1; ++n2) {
        // n2 * (n2 + 1)
        const float_type n2n2p1 = _an[n2 - 1];

        std::complex<float_type> this_J12, this_J21, this_RgJ12, this_RgJ21;
        // filter out half the components because of symmetry
        if (signs[n1 + n2 - 1] > 0) {
          for (uint_fast32_t ig = 1; ig < ngauss + 1; ++ig) {
            const float_type wigner_n1 = wigner_d(ig - 1, n1 - 1);
            const float_type dwigner_n1 = dwigner_d(ig - 1, n1 - 1);
            const float_type wigner_n2 = wigner_d(ig - 1, n2 - 1);
            const float_type dwigner_n2 = dwigner_d(ig - 1, n2 - 1);

            const float_type wn1dwn2 = wigner_n1 * dwigner_n2;
            const float_type dwn1wn2 = dwigner_n1 * wigner_n2;
            const float_type dwn1dwn2 = dwigner_n1 * dwigner_n2;

            const float_type jkrn1 = _jkr(ig - 1, n1 - 1);
            const float_type ykrn1 = _ykr(ig - 1, n1 - 1);
            // spherical Hankel function of the first kind
            const std::complex<float_type> hkrn1(jkrn1, ykrn1);
            const float_type djkrn1 = _djkr(ig - 1, n1 - 1);
            const float_type dykrn1 = _dykr(ig - 1, n1 - 1);
            const std::complex<float_type> dhkrn1(djkrn1, dykrn1);
            const std::complex<float_type> jkrmrn2 = _jkrmr(ig - 1, n2 - 1);
            const std::complex<float_type> djkrmrn2 = _djkrmr(ig - 1, n2 - 1);

            const std::complex<float_type> c1 = jkrmrn2 * jkrn1;
            const std::complex<float_type> b1 = jkrmrn2 * hkrn1;

            const std::complex<float_type> c2 = jkrmrn2 * djkrn1;
            const std::complex<float_type> b2 = jkrmrn2 * dhkrn1;

            const float_type krinvi = _krinv[ig - 1];
            const std::complex<float_type> c3 = krinvi * c1;
            const std::complex<float_type> b3 = krinvi * b1;

            const std::complex<float_type> c4 = jkrn1 * djkrmrn2;
            const std::complex<float_type> b4 = hkrn1 * djkrmrn2;

            const std::complex<float_type> krmrinvi = _krmrinv[ig - 1];
            const std::complex<float_type> c5 = c1 * krmrinvi;
            const std::complex<float_type> b5 = b1 * krmrinvi;

            const float_type wr2i = wr2[ig - 1];
            const float_type dr_over_ri = _dr_over_r[ig - 1];

            const float_type f1 = wr2i * dwn1dwn2;
            const float_type f2 = wr2i * dr_over_ri * n1n1p1 * wn1dwn2;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn2(krmr) *
            //    ([kr*hn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * hn1(kr) / kr
            this_J12 += f1 * b2 + f2 * b3;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn2(krmr) *
            //    ([kr*jn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * jn1(kr) / kr
            this_RgJ12 += f1 * c2 + f2 * c3;

            const float_type f3 = wr2i * dr_over_ri * n2n2p1 * dwn1wn2;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * hn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * hn1(kr) / krmr
            this_J21 += f1 * b4 + f3 * b5;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * jn1(kr) / krmr
            this_RgJ21 += f1 * c4 + f3 * c5;
          }
          // prefactor sqrt{(2n1+1)*(2n2+1)/[n1*(n1+1)*n2*(n2+1)]}
          const float_type an12 = 2. * _ann(n1 - 1, n2 - 1);
          J12(n1 - 1, n2 - 1) = an12 * this_J12;
          J21(n1 - 1, n2 - 1) = an12 * this_J21;
          RgJ12(n1 - 1, n2 - 1) = an12 * this_RgJ12;
          RgJ21(n1 - 1, n2 - 1) = an12 * this_RgJ21;
        }
      }
    }
    for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
      const uint_fast32_t k1 = n1;
      const uint_fast32_t kk1 = k1 + nmax;
      for (uint_fast32_t n2 = 1; n2 < nmax + 1; ++n2) {
        const uint_fast32_t k2 = n2;
        const uint_fast32_t kk2 = k2 + nmax;

        const std::complex<float_type> icompl(0., 1.);
        // no idea why we multiply with i: completely unnecessary...
        // (code also works if you leave out the i factor)
        // sign differences are due to a sign difference between the
        // implementation and documentation
        const std::complex<float_type> this_J12 = -icompl * J12(n1 - 1, n2 - 1);
        const std::complex<float_type> this_RgJ12 =
            -icompl * RgJ12(n1 - 1, n2 - 1);
        const std::complex<float_type> this_J21 = icompl * J21(n1 - 1, n2 - 1);
        const std::complex<float_type> this_RgJ21 =
            icompl * RgJ21(n1 - 1, n2 - 1);

        Q(k1 - 1, k2 - 1) = _k2mr * this_J21 + _k2 * this_J12;
        RgQ(k1 - 1, k2 - 1) = _k2mr * this_RgJ21 + _k2 * this_RgJ12;

        Q(kk1 - 1, kk2 - 1) = _k2mr * this_J12 + _k2 * this_J21;
        RgQ(kk1 - 1, kk2 - 1) = _k2mr * this_RgJ12 + _k2 * this_RgJ21;
      }
    }

    // func_TT
    Q.plu_inverse();

    for (uint_fast32_t i = 0; i < nmax; ++i) {
      const uint_fast32_t li = (i + 1) * (i + 2) - 1;
      for (uint_fast32_t j = 0; j < nmax; ++j) {
        const uint_fast32_t lj = (j + 1) * (j + 2) - 1;
        for (uint_fast32_t k = 0; k < nmax2; ++k) {
          _T(li, lj) -= RgQ(i, k) * Q(k, j);
          _T(_Lmax + li, lj) -= RgQ(nmax + i, k) * Q(k, j);
          _T(li, _Lmax + lj) -= RgQ(i, k) * Q(k, nmax + j);
          _T(_Lmax + li, _Lmax + lj) -= RgQ(nmax + i, k) * Q(k, nmax + j);
        }
      }
    }
  }

  /**
   * @brief Compute the missing elements of the T-matrix.
   */
  inline void compute_additional_elements() {

    for (uint_fast32_t m = 1; m < _nmax; ++m) {

      const float_type m2 = m * m;
      const uint_fast32_t nmax2 = 2 * _nmax;
      const uint_fast32_t nm = _nmax + 1 - m;
      const uint_fast32_t nm2 = 2 * nm;

      std::vector<int_fast8_t> signs(nmax2);
      Matrix<float_type> wigner_d(2 * _ngauss, _nmax);
      Matrix<float_type> dwigner_d(2 * _ngauss, _nmax);
      std::vector<float_type> wr2(_ngauss);
      Matrix<std::complex<float_type>> J11(_nmax, _nmax);
      Matrix<std::complex<float_type>> J12(_nmax, _nmax);
      Matrix<std::complex<float_type>> J21(_nmax, _nmax);
      Matrix<std::complex<float_type>> J22(_nmax, _nmax);
      Matrix<std::complex<float_type>> RgJ11(_nmax, _nmax);
      Matrix<std::complex<float_type>> RgJ12(_nmax, _nmax);
      Matrix<std::complex<float_type>> RgJ21(_nmax, _nmax);
      Matrix<std::complex<float_type>> RgJ22(_nmax, _nmax);
      Matrix<std::complex<float_type>> Q(nm2, nm2);
      Matrix<std::complex<float_type>> RgQ(nm2, nm2);
      std::vector<float_type> ds(_ngauss);
      std::vector<float_type> dss(_ngauss);

      int_fast8_t si = 1;
      for (uint_fast32_t n = 0; n < nmax2; ++n) {
        si = -si;
        signs[n] = si;
      }
      for (uint_fast32_t ig = 1; ig < _ngauss + 1; ++ig) {
        const uint_fast32_t i1 = _ngauss + ig;
        const uint_fast32_t i2 = _ngauss + 1 - ig;
        std::vector<float_type> dv1(_nmax), dv2(_nmax);
        SpecialFunctions::wigner_dn_0m(_costheta[i1 - 1], _nmax, m, &dv1[0],
                                       &dv2[0]);
        for (uint_fast32_t n = 0; n < _nmax; ++n) {
          si = signs[n];
          wigner_d(i1 - 1, n) = dv1[n];
          wigner_d(i2 - 1, n) = si * dv1[n];
          dwigner_d(i1 - 1, n) = dv2[n];
          dwigner_d(i2 - 1, n) = -si * dv2[n];
        }
      }
      for (uint_fast32_t ig = 0; ig < _ngauss; ++ig) {
        // move to class member
        wr2[ig] = _weights[ig] * _r2[ig];
        ds[ig] = _sinthetainv[ig] * m * wr2[ig];
        dss[ig] = _sintheta2inv[ig] * m2;
      }
      for (uint_fast32_t n1 = m; n1 < _nmax + 1; ++n1) {
        // n1 * (n1 + 1)
        const float_type n1n1p1 = _an[n1 - 1];
        for (uint_fast32_t n2 = m; n2 < _nmax + 1; ++n2) {
          // n2 * (n2 + 1)
          const float_type n2n2p1 = _an[n2 - 1];

          std::complex<float_type> this_J11, this_J12, this_J21, this_J22,
              this_RgJ11, this_RgJ12, this_RgJ21, this_RgJ22;
          const int_fast8_t si = signs[n1 + n2 - 1];
          for (uint_fast32_t ig = 0; ig < _ngauss; ++ig) {
            const float_type wigner_n1 = wigner_d(ig, n1 - 1);
            const float_type dwigner_n1 = dwigner_d(ig, n1 - 1);
            const float_type wigner_n2 = wigner_d(ig, n2 - 1);
            const float_type dwigner_n2 = dwigner_d(ig, n2 - 1);

            const float_type wn1wn2 = wigner_n1 * wigner_n2;
            const float_type wn1dwn2 = wigner_n1 * dwigner_n2;
            const float_type dwn1wn2 = dwigner_n1 * wigner_n2;

            const float_type jkrn1 = _jkr(ig, n1 - 1);
            const float_type ykrn1 = _ykr(ig, n1 - 1);
            // spherical Hankel function of the first kind
            const std::complex<float_type> hkrn1(jkrn1, ykrn1);
            const float_type djkrn1 = _djkr(ig, n1 - 1);
            const float_type dykrn1 = _dykr(ig, n1 - 1);
            const std::complex<float_type> dhkrn1(djkrn1, dykrn1);
            const std::complex<float_type> jkrmrn2 = _jkrmr(ig, n2 - 1);
            const std::complex<float_type> djkrmrn2 = _djkrmr(ig, n2 - 1);

            const std::complex<float_type> c1 = jkrmrn2 * jkrn1;
            const std::complex<float_type> b1 = jkrmrn2 * hkrn1;

            const std::complex<float_type> c2 = jkrmrn2 * djkrn1;
            const std::complex<float_type> b2 = jkrmrn2 * dhkrn1;

            const float_type krinvi = _krinv[ig];

            const std::complex<float_type> c4 = djkrmrn2 * jkrn1;
            const std::complex<float_type> b4 = djkrmrn2 * hkrn1;

            const std::complex<float_type> krmrinvi = _krmrinv[ig];

            const float_type dr_over_ri = _dr_over_r[ig];

            if (si < 0) {
              const float_type dsi = ds[ig];

              const std::complex<float_type> c6 = djkrmrn2 * djkrn1;
              const std::complex<float_type> b6 = djkrmrn2 * dhkrn1;

              const std::complex<float_type> c7 = c4 * krinvi;
              const std::complex<float_type> b7 = b4 * krinvi;

              const std::complex<float_type> c8 = c2 * krmrinvi;
              const std::complex<float_type> b8 = b2 * krmrinvi;

              const float_type e1 = dsi * (wn1dwn2 + dwn1wn2);
              // (m / sintheta) * jn2(krmr) * hn1(kr) *
              // (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
              this_J11 += e1 * b1;
              // (m / sintheta) * jn2(krmr) * jn1(kr) *
              // (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
              this_RgJ11 += e1 * c1;

              const float_type factor = dsi * dr_over_ri * wn1wn2;
              const float_type e2 = factor * n1n1p1;
              const float_type e3 = factor * n2n2p1;
              // (m / sintheta) * ([krmr*kn2(krmr)]'/krmr) * ([kr*hn1(kr)]'/kr)
              //    * (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
              // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n1 * (n1 + 1)
              //    * ([krmr*jn2(krmr)]'/krmr) * hn1(kr) / kr
              // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n2 * (n2 + 1)
              //    * jn2(krmr) * ([kr*hn1(kr]'/kr) / krmr
              this_J22 += e1 * b6 + e2 * b7 + e3 * b8;
              // (m / sintheta) * ([krmr*kn2(krmr)]'/krmr) * ([kr*jn1(kr)]'/kr)
              //    * (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
              // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n1 * (n1 + 1)
              //    * ([krmr*jn2(krmr)]'/krmr) * jn1(kr) / kr
              // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n2 * (n2 + 1)
              //    * jn2(krmr) * ([kr*jn1(kr]'/kr) / krmr
              this_RgJ22 += e1 * c6 + e2 * c7 + e3 * c8;
            } else {
              const float_type wr2i = wr2[ig];
              const std::complex<float_type> c3 = krinvi * c1;
              const std::complex<float_type> b3 = krinvi * b1;
              const std::complex<float_type> c5 = c1 * krmrinvi;
              const std::complex<float_type> b5 = b1 * krmrinvi;

              const float_type f1 =
                  wr2i * (wn1wn2 * dss[ig] + dwigner_n1 * dwigner_n2);
              const float_type f2 = wr2i * dr_over_ri * n1n1p1 * wn1dwn2;
              // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
              //        ddn1_0m/dtheta * ddn2_0m/dtheta) * jn2(krmr) *
              //    ([kr*hn1(kr)]'/kr)
              // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
              //    jn2(krmr) * hn1(kr) / kr
              this_J12 += f1 * b2 + f2 * b3;
              // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
              //        ddn1_0m/dtheta * ddn2_0m/dtheta) * jn2(krmr) *
              //    ([kr*jn1(kr)]'/kr)
              // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
              //    jn2(krmr) * jn1(kr) / kr
              this_RgJ12 += f1 * c2 + f2 * c3;

              const float_type f3 = wr2i * dr_over_ri * n2n2p1 * dwn1wn2;
              // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
              //        ddn1_0m/dtheta * ddn2_0m/dtheta) * hn1(kr) *
              //    ([krmr*jn2(krmr)]'/krmr)
              // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
              //    jn2(krmr) * hn1(kr) / krmr
              this_J21 += f1 * b4 + f3 * b5;
              // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
              //        ddn1_0m/dtheta * ddn2_0m/dtheta) * jn1(kr) *
              //    ([krmr*jn2(krmr)]'/krmr)
              // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
              //    jn2(krmr) * jn1(kr) / krmr
              this_RgJ21 += f1 * c4 + f3 * c5;
            }
          }
          // prefactor sqrt{(2n1+1)*(2n2+1)/[n1*(n1+1)*n2*(n2+1)]}
          const float_type an12 = 2. * _ann(n1 - 1, n2 - 1);
          J11(n1 - 1, n2 - 1) = this_J11 * an12;
          J12(n1 - 1, n2 - 1) = this_J12 * an12;
          J21(n1 - 1, n2 - 1) = this_J21 * an12;
          J22(n1 - 1, n2 - 1) = this_J22 * an12;
          RgJ11(n1 - 1, n2 - 1) = this_RgJ11 * an12;
          RgJ12(n1 - 1, n2 - 1) = this_RgJ12 * an12;
          RgJ21(n1 - 1, n2 - 1) = this_RgJ21 * an12;
          RgJ22(n1 - 1, n2 - 1) = this_RgJ22 * an12;
        }
      }
      for (uint_fast32_t n1 = m; n1 < _nmax + 1; ++n1) {
        const uint_fast32_t k1 = n1 + 1 - m;
        const uint_fast32_t kk1 = k1 + nm;
        for (uint_fast32_t n2 = m; n2 < _nmax + 1; ++n2) {
          const uint_fast32_t k2 = n2 + 1 - m;
          const uint_fast32_t kk2 = k2 + nm;

          const std::complex<float_type> icompl(0., 1.);
          // a factor -i is missing in J11 and J22
          // to compensate for this, we multiply J12 and J21 with -i too
          // we then multiply J11 and J22 with -1, so that it is wrong again?
          // not sure how to make sense of this...
          const std::complex<float_type> this_J11 = -J11(n1 - 1, n2 - 1);
          const std::complex<float_type> this_RgJ11 = -RgJ11(n1 - 1, n2 - 1);
          const std::complex<float_type> this_J12 =
              -icompl * J12(n1 - 1, n2 - 1);
          const std::complex<float_type> this_RgJ12 =
              -icompl * RgJ12(n1 - 1, n2 - 1);
          const std::complex<float_type> this_J21 =
              icompl * J21(n1 - 1, n2 - 1);
          const std::complex<float_type> this_RgJ21 =
              icompl * RgJ21(n1 - 1, n2 - 1);
          const std::complex<float_type> this_J22 = -J22(n1 - 1, n2 - 1);
          const std::complex<float_type> this_RgJ22 = -RgJ22(n1 - 1, n2 - 1);

          Q(k1 - 1, k2 - 1) = _k2mr * this_J21 + _k2 * this_J12;
          RgQ(k1 - 1, k2 - 1) = _k2mr * this_RgJ21 + _k2 * this_RgJ12;

          Q(k1 - 1, kk2 - 1) = _k2mr * this_J11 + _k2 * this_J22;
          RgQ(k1 - 1, kk2 - 1) = _k2mr * this_RgJ11 + _k2 * this_RgJ22;

          Q(kk1 - 1, k2 - 1) = _k2mr * this_J22 + _k2 * this_J11;
          RgQ(kk1 - 1, k2 - 1) = _k2mr * this_RgJ22 + _k2 * this_RgJ11;

          Q(kk1 - 1, kk2 - 1) = _k2mr * this_J12 + _k2 * this_J21;
          RgQ(kk1 - 1, kk2 - 1) = _k2mr * this_RgJ12 + _k2 * this_RgJ21;
        }
      }

      Q.plu_inverse();

      for (uint_fast32_t i = 0; i < nm; ++i) {
        const uint_fast32_t lip = (m + i) * (m + i + 1) + m - 1;
        const uint_fast32_t lim = (m + i) * (m + i + 1) - m - 1;
        for (uint_fast32_t j = 0; j < nm; ++j) {
          const uint_fast32_t ljp = (m + j) * (m + j + 1) + m - 1;
          const uint_fast32_t ljm = (m + j) * (m + j + 1) - m - 1;
          for (uint_fast32_t k = 0; k < nm2; ++k) {
            _T(lip, ljp) -= RgQ(i, k) * Q(k, j);
            _T(_Lmax + lip, ljp) -= RgQ(nm + i, k) * Q(k, j);
            _T(lip, _Lmax + ljp) -= RgQ(i, k) * Q(k, nm + j);
            _T(_Lmax + lip, _Lmax + ljp) -= RgQ(nm + i, k) * Q(k, nm + j);

            _T(lim, ljm) -= RgQ(i, k) * Q(k, j);
            _T(_Lmax + lim, ljm) -= RgQ(nm + i, k) * Q(k, j);
            _T(lim, _Lmax + ljm) -= RgQ(i, k) * Q(k, nm + j);
            _T(_Lmax + lim, _Lmax + ljm) -= RgQ(nm + i, k) * Q(k, nm + j);
          }
        }
      }
    }
  }

  /**
   * @brief Get a specific element of the T-matrix.
   *
   * The T-matrix consists of four blocks:
   * @f[
   *    T = \begin{pmatrix}
   *      T^{(11)} & T^{(12)} \\
   *      T^{(21)} & T^{(22)}
   *    \end{pmatrix}.
   * @f]
   * Each of these blocks is an @f$L_{max}\times{}L_{max}@f$ matrix indexed
   * using the combined index
   * @f[
   *    l = n (n+1) + m,
   * @f]
   * where @f$n = 1...n_{max}@f$ and @f$m=-n...n@f$. This notation is based on
   * Tsang, Kong & Shin, 1984, Radio Science, 19, 629
   * (https://doi.org/10.1029/RS019i002p00629).
   *
   * This function returns the element @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   *
   * @param i1 Row index of the desired T-matrix quarter, @f$i_1@f$.
   * @param n1 Order of the row index of the element, @f$n_1@f$.
   * @param m1 Degree of the row index of the element, @f$m_1@f$.
   * @param i2 Column index of the desired T-matrix quarter, @f$i_2@f$.
   * @param n2 Order of the column index of the element, @f$n_2@f$.
   * @param m2 Degree of the column index of the element, @f$m_2@f$.
   * @return Corresponding element, @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   */
  inline const std::complex<float_type> &
  operator()(const uint_fast8_t i1, const uint_fast32_t n1,
             const int_fast32_t m1, const uint_fast8_t i2,
             const uint_fast32_t n2, const int_fast32_t m2) const {

    const uint_fast32_t l1 = i1 * _Lmax + n1 * (n1 + 1) + m1 - 1;
    const uint_fast32_t l2 = i2 * _Lmax + n2 * (n2 + 1) + m2 - 1;
    return _T(l1, l2);
  }

  /**
   * @brief Get the scattering matrix for a scattering from the given input
   * angles to the given output angles at a particle with the given orientation.
   *
   * The scattering matrix @f$Z@f$ is a @f$4\times{}4@f$ matrix that transforms
   * the 4 elements of an incoming Stokes vector into the 4 elements of the
   * scattered Stokes vector after a scattering event. This matrix depends
   * on the direction of the incoming radiation, @f$\vec{n}_i@f$, the direction
   * of the scattered radiation, @f$\vec{n}_s@f$ (both measured in some
   * arbitrary fixed reference *laboratory* frame with its origin at the centre
   * of the particle), and the orientation of the scattering particle, generally
   * encoded in the three Euler angles that link the laboratory frame with a
   * frame fixed to the particle:
   *  - the angle @f$\alpha{}@f$ between the original @f$y_L@f$ axis and the
   *    projection of the new @f$y_P@f$ axis in the original @f$x_Ly_L@f$ plane,
   *    after a counterclockwise rotation along the @f$z_L@f$ axis.
   *  - the angle @f$\beta{}@f$ between the @f$z_L@f$ axis and the @f$z_P@f$
   *    axis after a counterclockwise rotation along the projection of the
   *    @f$y_P@f$ axis in the @f$x_Ly_L@f$ plane.
   *  - the angle @f$\gamma{}@f$ between the projection of the @f$y_P@f$ axis in
   *    the @f$x_Ly_L@f$ plane and the @f$y_P@f$ axis after a counterclockwise
   *    rotation along the @f$z_P@f$ axis.
   *
   * Because of the spheroidal symmetry of the scattering particle, we can
   * always choose the particle reference frame so that the @f$y_P@f$ axis and
   * its projection in the @f$x_Ly_L@f$ plane coincide, i.e. @f$\gamma{}=0@f$.
   * In this case, the angles @f$\alpha{}@f$ and @f$\beta{}@f$ will simply be
   * the azimuth and zenith angles for the @f$z_P@f$ axis, which we choose to
   * align with the rotation axis of the spheroidal particle. Similarly, we
   * can express the incoming and scattered directions in terms of a zenith
   * angle @f$\theta{}@f$ and an azimuth angle @f$\phi{}@f$. These are the 6
   * input parameters for the scattering function.
   *
   * The scattering matrix @f$Z@f$ is related to the @f$2\times{}2@f$
   * scattering matrix @f$S@f$ that transforms the two (complex) amplitudes
   * of the electromagnetic field perpendicular to the incoming direction into
   * the two amplitudes perpendicular to the outgoing direction. The relation
   * between the matrix @f$S@f$ and the matrix @f$Z@f$ is given in section VI
   * of chapter 1 of Mishchenko, Hovenier & Travis, 2000, Light Scattering by
   * Nonspherical Particles
   * (https://www.elsevier.com/books/light-scattering-by-nonspherical-particles/mishchenko/978-0-12-498660-2),
   * where @f$X_{ij}@f$ represents the element on row @f$i@f$, column @f$j@f$ of
   * the matrix @f$X@f$:
   * @f[
   *    Z_{11} = \frac{1}{2} \left(
   *      |S_{11}|^2 + |S_{12}|^2 + |S_{21}|^2 + |S_{22}|^2
   *    \right),
   * @f]
   * @f[
   *    Z_{12} = \frac{1}{2} \left(
   *      |S_{11}|^2 - |S_{12}|^2 + |S_{21}|^2 - |S_{22}|^2
   *    \right),
   * @f]
   * @f[
   *    Z_{13} = -\Re\left( S_{11} S^*_{12} + S_{22} S^*_{21} \right),
   * @f]
   * @f[
   *    Z_{14} = -\Im\left( S_{11} S^*_{12} - S_{22} S^*_{21} \right),
   * @f]
   * @f[
   *    Z_{21} = \frac{1}{2} \left(
   *      |S_{11}|^2 + |S_{12}|^2 - |S_{21}|^2 - |S_{22}|^2
   *    \right),
   * @f]
   * @f[
   *    Z_{22} = \frac{1}{2} \left(
   *      |S_{11}|^2 - |S_{12}|^2 - |S_{21}|^2 + |S_{22}|^2
   *    \right),
   * @f]
   * @f[
   *    Z_{23} = -\Re\left( S_{11} S^*_{12} - S_{22} S^*_{21} \right),
   * @f]
   * @f[
   *    Z_{24} = -\Im\left( S_{11} S^*_{12} + S_{22} S^*_{21} \right),
   * @f]
   * @f[
   *    Z_{31} = -\Re\left( S_{11} S^*_{21} + S_{22} S^*_{12} \right),
   * @f]
   * @f[
   *    Z_{32} = -\Re\left( S_{11} S^*_{21} - S_{22} S^*_{12} \right),
   * @f]
   * @f[
   *    Z_{33} = \Re\left( S_{11} S^*_{22} + S_{12} S^*_{21} \right),
   * @f]
   * @f[
   *    Z_{34} = \Im\left( S_{11} S^*_{22} + S_{21} S^*_{12} \right),
   * @f]
   * @f[
   *    Z_{41} = -\Im\left( S_{21} S^*_{11} + S_{22} S^*_{12} \right),
   * @f]
   * @f[
   *    Z_{42} = -\Im\left( S_{21} S^*_{11} - S_{22} S^*_{12} \right),
   * @f]
   * @f[
   *    Z_{43} = \Im\left( S_{22} S^*_{11} - S_{12} S^*_{21} \right),
   * @f]
   * @f[
   *    Z_{44} = \Re\left( S_{22} S^*_{11} - S_{12} S^*_{21} \right),
   * @f]
   * where @f$\Re(z)@f$ and @f$\Im(z)@f$ represent respectively the real part
   * @f$x@f$ and imaginary part @f$y@f$ of the complex number @f$z = x + iy@f$.
   *
   * The scattering matrix @f$S@f$ can be computed in the particle reference
   * frame from the T matrix (Mishchenko, 2000, Applied Optics, 39, 1026;
   * https://doi.org/10.1364/AO.39.001026):
   * @f[
   *    S^P_{11} = \frac{1}{k} \sum_{n=1}^{n_{max}} \sum_{n'=1}^{n_{max}}
   *      \sum_{m = -\min(n, n')}^{\min(n, n')}
   *      i^{n' - n -1} \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
   *      e^{im(\phi{}_s^P - \phi{}_i^P)} \\\left(
   *        T^{11}_{mnmn'} \frac{m^2}{\sin^2(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{21}_{mnmn'} \frac{m}{\sin(\theta{})}
   *          \frac{d}{d\theta{}} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{12}_{mnmn'} \frac{m}{\sin(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{22}_{mnmn'} \frac{d}{d\theta{}} d^{n}_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P)
   *      \right),
   * @f]
   * @f[
   *    S^P_{12} = -\frac{i}{k} \sum_{n=1}^{n_{max}} \sum_{n'=1}^{n_{max}}
   *      \sum_{m = -\min(n, n')}^{\min(n, n')}
   *      i^{n' - n -1} \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
   *      e^{im(\phi{}_s^P - \phi{}_i^P)} \\\left(
   *        T^{11}_{mnmn'} \frac{m}{\sin(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{21}_{mnmn'} \frac{d}{d\theta{}} d^{n}_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{12}_{mnmn'} \frac{m^2}{\sin^2(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{22}_{mnmn'} \frac{m}{\sin(\theta{})}
   *          \frac{d}{d\theta{}} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P)
   *      \right),
   * @f]
   * @f[
   *    S^P_{21} = \frac{i}{k} \sum_{n=1}^{n_{max}} \sum_{n'=1}^{n_{max}}
   *      \sum_{m = -\min(n, n')}^{\min(n, n')}
   *      i^{n' - n -1} \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
   *      e^{im(\phi{}_s^P - \phi{}_i^P)} \\\left(
   *        T^{11}_{mnmn'} \frac{m}{\sin(\theta{})}
   *          \frac{d}{d\theta{}} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{21}_{mnmn'} \frac{m^2}{\sin^2(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{12}_{mnmn'} \frac{d}{d\theta{}} d^{n}_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{22}_{mnmn'} \frac{m}{\sin(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P)
   *      \right),
   * @f]
   * @f[
   *    S^P_{22} = \frac{1}{k} \sum_{n=1}^{n_{max}} \sum_{n'=1}^{n_{max}}
   *      \sum_{m = -\min(n, n')}^{\min(n, n')}
   *      i^{n' - n -1} \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}
   *      e^{im(\phi{}_s^P - \phi{}_i^P)} \\\left(
   *        T^{11}_{mnmn'} \frac{d}{d\theta{}} d^{n}_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{21}_{mnmn'} \frac{m}{\sin(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          \frac{d}{d\theta{}} d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{12}_{mnmn'} \frac{m}{\sin(\theta{})}
   *          \frac{d}{d\theta{}} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P) +
   *        T^{22}_{mnmn'} \frac{m^2}{\sin^2(\theta{})} d^n_{0m}(\theta{}_s^P)
   *          d^{n'}_{0m}(\theta{}_i^P)
   *      \right).
   * @f]
   * In this expression, the angles @f$\phi{}_i^P, \phi{}_s^P, \theta{}_i^P@f$
   * and @f$\theta_s^P@f$ represent the azimuth and zenith angles of the
   * incoming and outgoing ray in the particle reference frame. Note that
   * the scattering matrix in this frame strictly speaking only depends on
   * three angles: the zenith angles @f$\theta{}_i^P@f$ and @f$\theta{}_s^P@f$
   * and the relative azimuth angle @f$\phi{}_s^P - \phi{}_i^P@f$.
   *
   * Also note that we use symbols similar to the ones introduced in Mishchenko
   * (2000) in the code:
   * @f[
   *    c_{nn'} = i^{n'-n-1} \sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)},
   * @f]
   * @f[
   *    \pi{}_{mn}(\theta{}) = \frac{m}{\sin(\theta{})} d^n_{0m}(\theta{}),
   * @f]
   * and
   * @f[
   *    \tau{}_{mn}(\theta{}) = \frac{d}{d\theta{}} d^n_{0m}(\theta{}).
   * @f]
   *
   * To compute the @f$Z@f$ matrix in the original laboratory reference frame,
   * we need to rotate the scattering matrix @f$S^P@f$ into the laboratory
   * frame to get @f$S^L@f$. This requires knowledge of the reference frame
   * transformation, and hence the rotations encoded by @f$\alpha{}@f$ and
   * @f$\beta{}@f$. If we denote the unit vectors along the Cartesian coordinate
   * axes in the laboratory frame as @f$\hat{x}_L, \hat{y}_L, \hat{z}_L@f$ and
   * the unit vectors in the particle frame as @f$\hat{x}_P, \hat{y}_P,
   * \hat{z}_P@f$, and we also introduce the intermediate frame with unit
   * vectors @f$\hat{x}_I, \hat{y}_I, \hat{z}_I@f$ that represents the frame
   * after the first rotation with angle @f$\alpha{}@f$, then we can derive
   * the following relations between the various unit vectors:
   * @f[
   *    \begin{cases}
   *      \hat{x}_I = \cos(\alpha{}) \hat{x}_L + \sin(\alpha{}) \hat{y}_L, \\
   *      \hat{y}_I = -\sin(\alpha{}) \hat{x}_L + \cos(\alpha{}) \hat{y}_L, \\
   *      \hat{z}_I = \hat{z}_L,
   *    \end{cases}
   * @f]
   * @f[
   *    \begin{cases}
   *      \hat{x}_L = \cos(\alpha{}) \hat{x}_I - \sin(\alpha{}) \hat{y}_I, \\
   *      \hat{y}_L = \sin(\alpha{}) \hat{x}_I + \cos(\alpha{}) \hat{y}_I, \\
   *      \hat{z}_L = \hat{z}_I,
   *    \end{cases}
   * @f]
   * @f[
   *    \begin{cases}
   *      \hat{x}_P = \cos(\beta{}) \hat{x}_I - \sin(\beta{}) \hat{z}_I =
   *        \cos(\alpha{}) \cos(\beta{}) \hat{x}_L +
   *        \sin(\alpha{}) \cos(\beta{}) \hat{y}_L - \sin(\beta{}) \hat{z}_L, \\
   *      \hat{y}_P = \hat{y}_I =
   *        -\sin(\alpha{}) \hat{x}_L + \cos(\alpha{}) \hat{y}_L, \\
   *      \hat{z}_P = \sin(\beta{}) \hat{x}_I + \cos(\beta{}) \hat{z}_I =
   *        \cos(\alpha{}) \sin(\beta{}) \hat{x}_L +
   *        \sin(\alpha{}) \sin(\beta{}) \hat{y}_L + \cos(\beta{}) \hat{z}_L,
   *    \end{cases}
   * @f]
   * @f[
   *    \begin{cases}
   *      \hat{x}_L = \cos(\alpha{}) \cos(\beta{}) \hat{x}_P
   *        - \sin(\alpha{}) \hat{y}_P
   *        + \cos(\alpha{}) \sin(\beta{}) \hat{z}_P, \\
   *      \hat{y}_L = \sin(\alpha{}) \cos(\beta{}) \hat{x}_P
   *        + \cos(\alpha{}) \hat{y}_P
   *        + \sin(\alpha{}) \sin(\beta{}) \hat{z}_P, \\
   *      \hat{z}_L = -\sin(\beta{}) \hat{x}_P + \cos(\beta{}) \hat{z}_P.
   *    \end{cases}
   * @f]
   *
   * The direction vectors @f$\vec{n}_i@f$ and @f$\vec{n}_s@f$ can be expanded
   * in both reference frames as
   * @f[
   *    \vec{n} = \sin(\theta{}) \cos(\phi{}) \hat{x}
   *      + \sin(\theta{}) \sin(\phi{}) \hat{y}
   *      + \cos(\theta{}) \hat{z},
   * @f]
   * which leads to the following transformation formulas for @f$\theta{}^P@f$
   * and @f$\phi{}^P@f$:
   * @f[
   *    \cos(\theta{}^P) = \cos(\theta{}^L) \cos(\beta{})
   *      + \sin(\theta{}^L) \sin(\beta{}) \left(
   *        \cos(\phi{}^L) \cos(\alpha{}) + \sin(\phi{}^L) \sin(\alpha{})
   *      \right) = \cos(\theta{}^L) \cos(\beta{})
   *      + \sin(\theta{}^L) \sin(\beta{}) \cos(\phi{}^L - \alpha{}),
   * @f]
   * @f[
   *    \cos(\phi{}^P) = \frac{1}{\sin(\theta{}^P)} \left[
   *      \sin(\theta{}^L) \cos(\beta{}) \left(
   *        \cos(\phi{}^L) \cos(\alpha{}) + \sin(\phi{}^L) \sin(\alpha{})
   *      \right) - \cos(\theta{}^L) \sin(\beta{})
   *    \right] =
   *    \frac{1}{\sin(\theta{}^P)} \left(
   *      \sin(\theta{}^L) \cos(\beta{}) \cos(\phi{}^L - \alpha{})
   *        - \cos(\theta{}^L) \sin(\beta{})
   *    \right),
   * @f]
   * @f[
   *    \sin(\phi{}^P) = \frac{1}{\sin(\theta{}^P)} \left[
   *      \sin(\theta{}^L) \left(
   *        \sin(\phi{}^L) \cos(\alpha{}) - \cos(\phi{}^L) \sin(\alpha{})
   *      \right)
   *    \right] = \frac{1}{\sin(\theta{}^P)} \left(
   *      \sin(\theta{}^L) \sin(\phi{}^L - \alpha{})
   *    \right).
   * @f]
   * These formulas could potentially cause problems if @f$\theta{}^P = 0@f$,
   * since then the division in the formulas for @f$\phi{}^P@f$ will go wrong.
   * The original code of Mishchenko does not explicitly address these issues.
   * Since @f$\theta{}^P = 0@f$ signals a degeneracy in @f$\phi{}^P@f$, we
   * are free to choose any arbitrary value for @f$\phi{}^P@f$ in this case.
   * We will conveniently choose @f$\phi{}^P=0@f$. Note that we never explicitly
   * require @f$\theta{}^P@f$ or @f$\phi{}^P@f$; we always require either their
   * @f$\sin{}@f$ or @f$\cos{}@f$ or expressions like
   * @f[
   *    \cos(\phi{}^P_s - \phi{}^P_i) =
   *      \cos(\phi{}^P_s) \cos(\phi{}^P_i) + \sin(\phi{}^P_s) \sin(\phi{}^P_i),
   * @f]
   * or
   * @f[
   *    \sin(\phi{}^P_s - \phi{}^P_i) =
   *      \sin(\phi{}^P_s) \cos(\phi{}^P_i) - \cos(\phi{}^P_s) \sin(\phi{}^P_i),
   * @f]
   * which are trivially calculated. To compute the azimuthal factor in the
   * scattering matrix expression, we recursively compute
   * @f[
   *    e^{im(\phi{}^P_s - \phi{}^P_i)} = \left(
   *      \cos(\phi{}^P_s - \phi{}^P_i) + i \sin(\phi{}^P_s - \phi{}^P_i)
   *    \right)^m.
   * @f]
   * This is much faster than having to first explicitly compute the angles
   * @f$\phi{}^P@f$ and then evaluating the trigonometric functions for every
   * value of @f$m@f$.
   *
   * @param alpha_radians Azimuth angle of the particle's rotation axis,
   * @f$\alpha{}@f$ (in radians).
   * @param beta_radians Zenith angle of the particle's rotation axis,
   * @f$\beta{}@f$ (in radians).
   * @param theta_in_radians Zenith angle of the incoming photon,
   * @f$\theta{}_i@f$ (in radians).
   * @param phi_in_radians Azimuth angle of the incoming photon, @f$\phi{}_i@f$
   * (in radians).
   * @param theta_out_radians Zenith angle of the scattered photon,
   * @f$\theta{}_s@f$ (in radians).
   * @param phi_out_radians Azimuth angle fo the scattered photon,
   * @f$\phi{}_s@f$ (in radians).
   * @return Scattering matrix for this scattering event.
   */
  inline Matrix<float_type> get_scattering_matrix(
      const float_type alpha_radians, const float_type beta_radians,
      const float_type theta_in_radians, const float_type phi_in_radians,
      const float_type theta_out_radians,
      const float_type phi_out_radians) const {

    // Mishchenko includes some (buggy) corrections for small angles
    // might be worth looking into this in a later stage...

    // compute all sines and cosines in one go; we need all of them anyway
    const float_type cosalpha = cos(alpha_radians);
    const float_type sinalpha = sin(alpha_radians);
    const float_type cosbeta = cos(beta_radians);
    const float_type sinbeta = sin(beta_radians);
    const float_type costheta_l_in = cos(theta_in_radians);
    const float_type sintheta_l_in = sin(theta_in_radians);
    const float_type costheta_l_out = cos(theta_out_radians);
    const float_type sintheta_l_out = sin(theta_out_radians);
    const float_type cosphi_l_in = cos(phi_in_radians);
    const float_type sinphi_l_in = sin(phi_in_radians);
    const float_type cosphi_l_out = cos(phi_out_radians);
    const float_type sinphi_l_out = sin(phi_out_radians);

    const float_type cosphirel_in =
        cosphi_l_in * cosalpha + sinphi_l_in * sinalpha;
    const float_type sinphirel_in =
        sinphi_l_in * cosalpha - cosphi_l_in * sinalpha;

    const float_type costheta_p_in =
        costheta_l_in * cosbeta + sintheta_l_in * sinbeta * cosphirel_in;
    const float_type sintheta_p_in =
        sqrt((1. - costheta_p_in) * (1. + costheta_p_in));
    float_type cosphi_p_in, sinphi_p_in, sintheta_p_in_inv;
    if (sintheta_p_in != 0.) {
      sintheta_p_in_inv = 1. / sintheta_p_in;
      cosphi_p_in =
          (cosbeta * sintheta_l_in * cosphirel_in - sinbeta * costheta_l_in) *
          sintheta_p_in_inv;
      sinphi_p_in = sintheta_l_in * sinphirel_in * sintheta_p_in_inv;
    } else {
      cosphi_p_in = 1.;
      sinphi_p_in = 0.;
      // this value will not be used, but we set it to something anyway
      sintheta_p_in_inv = 9000.;
    }

    const float_type cosphirel_out =
        cosphi_l_out * cosalpha + sinphi_l_out * sinalpha;
    const float_type sinphirel_out =
        sinphi_l_out * cosalpha - cosphi_l_out * sinalpha;

    const float_type costheta_p_out =
        costheta_l_out * cosbeta + sintheta_l_out * sinbeta * cosphirel_out;
    const float_type sintheta_p_out =
        sqrt((1. - costheta_p_out) * (1. + costheta_p_out));
    float_type cosphi_p_out, sinphi_p_out, sintheta_p_out_inv;
    if (sintheta_p_out != 0.) {
      sintheta_p_out_inv = 1. / sintheta_p_out;
      cosphi_p_out = (cosbeta * sintheta_l_out * cosphirel_out -
                      sinbeta * costheta_l_out) *
                     sintheta_p_out_inv;
      sinphi_p_out = sintheta_l_out * sinphirel_out * sintheta_p_out_inv;
    } else {
      cosphi_p_out = 1.;
      sinphi_p_out = 0.;
      // this value will not be used, but we set it to something anyway
      sintheta_p_out_inv = 9000.;
    }

    Matrix<float_type> B(3, 3);
    B(0, 0) = cosalpha * cosbeta;
    B(0, 1) = sinalpha * cosbeta;
    B(0, 2) = -sinbeta;
    B(1, 0) = -sinalpha;
    B(1, 1) = cosalpha;
    // B(1,2) remains 0.
    B(2, 0) = cosalpha * sinbeta;
    B(2, 1) = sinalpha * sinbeta;
    B(2, 2) = cosbeta;

    Matrix<float_type> AL_in(3, 2);
    AL_in(0, 0) = costheta_l_in * cosphi_l_in;
    AL_in(0, 1) = -sinphi_l_in;
    AL_in(1, 0) = costheta_l_in * sinphi_l_in;
    AL_in(1, 1) = cosphi_l_in;
    AL_in(2, 0) = -sintheta_l_in;
    // AL_in(2,1) remains 0.

    Matrix<float_type> AL_out(3, 2);
    AL_out(0, 0) = costheta_l_out * cosphi_l_out;
    AL_out(0, 1) = -sinphi_l_out;
    AL_out(1, 0) = costheta_l_out * sinphi_l_out;
    AL_out(1, 1) = cosphi_l_out;
    AL_out(2, 0) = -sintheta_l_out;
    // AL_out(2,1) remains 0.

    Matrix<float_type> AP_in(2, 3);
    AP_in(0, 0) = costheta_p_in * cosphi_p_in;
    AP_in(0, 1) = costheta_p_in * sinphi_p_in;
    AP_in(0, 2) = -sintheta_p_in;
    AP_in(1, 0) = -sinphi_p_in;
    AP_in(1, 1) = cosphi_p_in;
    // AP_in(1,2) remains 0.

    Matrix<float_type> AP_out(2, 3);
    AP_out(0, 0) = costheta_p_out * cosphi_p_out;
    AP_out(0, 1) = costheta_p_out * sinphi_p_out;
    AP_out(0, 2) = -sintheta_p_out;
    AP_out(1, 0) = -sinphi_p_out;
    AP_out(1, 1) = cosphi_p_out;
    // AP_out(1,2) remains 0.

    Matrix<float_type> C(3, 2);
    for (uint_fast8_t i = 0; i < 3; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          C(i, j) += B(i, k) * AL_in(k, j);
        }
      }
    }
    Matrix<float_type> R_in(2, 2);
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          R_in(i, j) += AP_in(i, k) * C(k, j);
        }
      }
    }

    for (uint_fast8_t i = 0; i < 3; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        C(i, j) = 0.;
        for (uint_fast8_t k = 0; k < 3; ++k) {
          C(i, j) += B(i, k) * AL_out(k, j);
        }
      }
    }
    Matrix<float_type> R_out(2, 2);
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          R_out(i, j) += AP_out(i, k) * C(k, j);
        }
      }
    }

    const float_type d =
        1. / (R_out(0, 0) * R_out(1, 1) - R_out(0, 1) * R_out(1, 0));
    const float_type temp = R_out(0, 0);
    R_out(0, 0) = R_out(1, 1) * d;
    R_out(0, 1) = -R_out(0, 1) * d;
    R_out(1, 0) = -R_out(1, 0) * d;
    R_out(1, 1) = temp * d;

    const std::complex<float_type> icompl(0., 1.);
    Matrix<std::complex<float_type>> c(_nmax, _nmax);
    std::complex<float_type> icomp_pow_nn = icompl;
    for (uint_fast32_t nn = 1; nn < _nmax + 1; ++nn) {
      std::complex<float_type> icomp_pow_m_n_m_1(-1.);
      for (uint_fast32_t n = 1; n < _nmax + 1; ++n) {
        // icomp_pow_nn*icomp_pow_m_n_m_1 now equals i^(nn - n - 1)
        c(n - 1, nn - 1) = icomp_pow_m_n_m_1 * icomp_pow_nn *
                           float_type(sqrt((2. * n + 1.) * (2. * nn + 1.) /
                                           (n * nn * (n + 1.) * (nn + 1.))));
        icomp_pow_m_n_m_1 /= icompl;
      }
      icomp_pow_nn *= icompl;
    }

    const std::complex<float_type> expiphi_p_out_m_in(
        cosphi_p_out * cosphi_p_in + sinphi_p_out * sinphi_p_in,
        sinphi_p_out * cosphi_p_in - cosphi_p_out * sinphi_p_in);
    std::complex<float_type> expimphi_p_out_m_in(1., 0.);
    std::complex<float_type> S11, S12, S21, S22;
    const TMatrix &T = *this;
    for (uint_fast32_t m = 0; m < _nmax; ++m) {
      const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));

      std::vector<float_type> pi_in(_nmax), tau_in(_nmax);
      SpecialFunctions::wigner_dn_0m_sinx(costheta_p_in, sintheta_p_in,
                                          sintheta_p_in_inv, _nmax, m,
                                          &pi_in[0], &tau_in[0]);
      std::vector<float_type> pi_out(_nmax), tau_out(_nmax);
      SpecialFunctions::wigner_dn_0m_sinx(costheta_p_out, sintheta_p_out,
                                          sintheta_p_out_inv, _nmax, m,
                                          &pi_out[0], &tau_out[0]);

      // we get the real and imaginary part of e^{im\phi{}} and multiply with
      // 2 to account for both m and -m
      const float_type fcos = 2. * expimphi_p_out_m_in.real();
      const float_type fsin = 2. * expimphi_p_out_m_in.imag();
      // recurse the exponential for the next iteration
      expimphi_p_out_m_in *= expiphi_p_out_m_in;

      for (uint_fast32_t nn = nmin; nn < _nmax + 1; ++nn) {
        const float_type pi_nn = m * pi_in[nn - 1];
        const float_type tau_nn = tau_in[nn - 1];
        for (uint_fast32_t n = nmin; n < _nmax + 1; ++n) {

          const float_type pi_n = m * pi_out[n - 1];
          const float_type tau_n = tau_out[n - 1];

          const std::complex<float_type> c_nnn = c(n - 1, nn - 1);

          const std::complex<float_type> T11nmnnm = T(0, n, m, 0, nn, m);
          const std::complex<float_type> T22nmnnm = T(1, n, m, 1, nn, m);
          if (m == 0) {
            const std::complex<float_type> factor = c_nnn * tau_n * tau_nn;
            S11 += factor * T22nmnnm;
            S22 += factor * T11nmnnm;
          } else {
            const std::complex<float_type> T12nmnnm = T(0, n, m, 1, nn, m);
            const std::complex<float_type> T21nmnnm = T(1, n, m, 0, nn, m);

            const std::complex<float_type> real_factor = c_nnn * fcos;
            const std::complex<float_type> imag_factor = c_nnn * fsin;

            const float_type pi_pi = pi_n * pi_nn;
            const float_type pi_tau = pi_n * tau_nn;
            const float_type tau_pi = tau_n * pi_nn;
            const float_type tau_tau = tau_n * tau_nn;

            S11 += real_factor * (T11nmnnm * pi_pi + T21nmnnm * tau_pi +
                                  T12nmnnm * pi_tau + T22nmnnm * tau_tau);
            S12 += imag_factor * (T11nmnnm * pi_tau + T21nmnnm * tau_tau +
                                  T12nmnnm * pi_pi + T22nmnnm * tau_pi);
            S21 -= imag_factor * (T11nmnnm * tau_pi + T21nmnnm * pi_pi +
                                  T12nmnnm * tau_tau + T22nmnnm * pi_tau);
            S22 += real_factor * (T11nmnnm * tau_tau + T21nmnnm * pi_tau +
                                  T12nmnnm * tau_pi + T22nmnnm * pi_pi);
          }
        }
      }
    }
    const float_type kinv = 1. / _k;
    S11 *= kinv;
    S12 *= kinv;
    S21 *= kinv;
    S22 *= kinv;

    const std::complex<float_type> cS11 = S11 * R_in(0, 0) + S12 * R_in(1, 0);
    const std::complex<float_type> cS12 = S11 * R_in(0, 1) + S12 * R_in(1, 1);
    const std::complex<float_type> cS21 = S21 * R_in(0, 0) + S22 * R_in(1, 0);
    const std::complex<float_type> cS22 = S21 * R_in(0, 1) + S22 * R_in(1, 1);

    S11 = R_out(0, 0) * cS11 + R_out(0, 1) * cS21;
    S12 = R_out(0, 0) * cS12 + R_out(0, 1) * cS22;
    S21 = R_out(1, 0) * cS11 + R_out(1, 1) * cS21;
    S22 = R_out(1, 0) * cS12 + R_out(1, 1) * cS22;

    ctm_warning("S11: %g + i%g", double(S11.real()), double(S11.imag()));
    ctm_warning("S12: %g + i%g", double(S12.real()), double(S12.imag()));
    ctm_warning("S21: %g + i%g", double(S21.real()), double(S21.imag()));
    ctm_warning("S22: %g + i%g", double(S22.real()), double(S22.imag()));

    Matrix<float_type> Z(4, 4);

    const float_type half(0.5);
    Z(0, 0) = (half * (S11 * conj(S11) + S12 * conj(S12) + S21 * conj(S21) +
                       S22 * conj(S22)))
                  .real();
    Z(0, 1) = (half * (S11 * conj(S11) - S12 * conj(S12) + S21 * conj(S21) -
                       S22 * conj(S22)))
                  .real();
    Z(0, 2) = (-S11 * conj(S12) - S22 * conj(S21)).real();
    Z(0, 3) = (icompl * (S11 * conj(S12) - S22 * conj(S21))).real();

    Z(1, 0) = (half * (S11 * conj(S11) + S12 * conj(S12) - S21 * conj(S21) -
                       S22 * conj(S22)))
                  .real();
    Z(1, 1) = (half * (S11 * conj(S11) - S12 * conj(S12) - S21 * conj(S21) +
                       S22 * conj(S22)))
                  .real();
    Z(1, 2) = (-S11 * conj(S12) + S22 * conj(S21)).real();
    Z(1, 3) = (icompl * (S11 * conj(S12) + S22 * conj(S21))).real();

    Z(2, 0) = (-S11 * conj(S21) - S22 * conj(S12)).real();
    Z(2, 1) = (-S11 * conj(S21) + S22 * conj(S12)).real();
    Z(2, 2) = (S11 * conj(S22) + S12 * conj(S21)).real();
    Z(2, 3) = (-icompl * (S11 * conj(S22) + S21 * conj(S12))).real();

    Z(3, 0) = (icompl * (S21 * conj(S11) + S22 * conj(S12))).real();
    Z(3, 1) = (icompl * (S21 * conj(S11) - S22 * conj(S12))).real();
    Z(3, 2) = (-icompl * (S22 * conj(S11) - S12 * conj(S21))).real();
    Z(3, 3) = (S22 * conj(S11) - S12 * conj(S21)).real();

    return Z;
  }
};

#endif // TMATRIX_HPP

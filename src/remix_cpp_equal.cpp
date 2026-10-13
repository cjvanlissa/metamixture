// Finite mixtures of random-effects meta-analysis models.
// yi | component k ~ N(mu[k], vi + tau2[k]), with fixed, known vi > 0.
//
// Place this file in metamixture/src/remix.cpp and run Rcpp::compileAttributes().
// Only Rcpp is required. Joint refinement calls stats::optim (L-BFGS-B).
//
// Use the C++ standard selected by R; this file does not request C++11.
// Regenerate RcppExports with Rcpp::compileAttributes() after replacing it.
// The exported R function is remix_cpp(), with no seed argument.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace metamixture_core {

const double INF = std::numeric_limits<double>::infinity();
const double LOG2PI = 1.837877066409345483560659472811;

inline void interrupt() {
    Rcpp::checkUserInterrupt();
}
inline double clamp(double x, double a, double b) {
    return std::max(a, std::min(b, x));
}

struct Data {
    std::vector<double> y, v, vknots;
    double center, scale, ymin, ymax, upper, vmin;
    int n;
    Data(const std::vector<double>& yi, const std::vector<double>& vi) {
        n = static_cast<int>(yi.size());
        if (n < 2 || vi.size() != yi.size())
            throw std::invalid_argument("yi and vi must have equal length >= 2.");
        for (int i = 0; i < n; ++i)
            if (!std::isfinite(yi[i]) || !std::isfinite(vi[i]) || vi[i] <= 0)
                throw std::invalid_argument("yi must be finite; vi must be finite and > 0.");
        double lo = *std::min_element(yi.begin(), yi.end());
        double hi = *std::max_element(yi.begin(), yi.end());
        center = 0.5 * lo + 0.5 * hi;
        double dev = std::max(std::fabs(lo - center), std::fabs(hi - center));
        scale = std::max(dev, std::sqrt(*std::max_element(vi.begin(), vi.end())));
        if (!std::isfinite(scale) || scale <= 0)
            throw std::invalid_argument("Data cannot be scaled safely in double precision.");
        y.resize(n); v.resize(n);
        for (int i = 0; i < n; ++i) {
            y[i] = (yi[i] - center) / scale;
            v[i] = (vi[i] / scale) / scale;
            if (!std::isfinite(y[i]) || !std::isfinite(v[i]) || v[i] <= 0)
                throw std::invalid_argument("Dynamic range is too large for double precision.");
        }
        ymin = *std::min_element(y.begin(), y.end());
        ymax = *std::max_element(y.begin(), y.end());
        upper = (ymax - ymin) * (ymax - ymin);
        vmin = *std::min_element(v.begin(), v.end());
        std::vector<double> sv = v;
        std::sort(sv.begin(), sv.end());
        int nk = std::min(n, 25);
        for (int j = 0; j < nk; ++j) {
            double s = sv[static_cast<std::size_t>(j) * (n - 1) / (nk - 1)];
            if (s > 0 && s < upper) vknots.push_back(s);
        }
    }
    double original_ll(double ll) const { return ll - n * std::log(scale); }
};

struct Params {
    std::vector<double> mu, tau, pi;
    explicit Params(int K = 0) : mu(K), tau(K), pi(K) {}
    int size() const { return static_cast<int>(mu.size()); }
};

// Log-sum-exp likelihood. Optional log weights avoid underflow in the
// unconstrained logits used by the joint optimizer. r is stored row-major.
double estep(const Data& d, const Params& p, std::vector<double>& r,
             const std::vector<double>* supplied_logpi = nullptr) {
    const int K = p.size();
    r.assign(static_cast<std::size_t>(d.n) * K, 0.0);
    std::vector<double> lp(K), a(K);
    for (int k = 0; k < K; ++k)
        lp[k] = supplied_logpi ? (*supplied_logpi)[k] :
            (p.pi[k] > 0 ? std::log(p.pi[k]) : -INF);
    long double ll = 0;
    for (int i = 0; i < d.n; ++i) {
        double mx = -INF;
        for (int k = 0; k < K; ++k) {
            double s = d.v[i] + p.tau[k];
            long double z = static_cast<long double>(d.y[i]) - p.mu[k];
            a[k] = lp[k] - 0.5 * (LOG2PI + std::log(s) +
                       static_cast<double>(z * z / s));
            mx = std::max(mx, a[k]);
        }
        if (!std::isfinite(mx)) return -INF;
        double den = 0;
        for (int k = 0; k < K; ++k) den += std::exp(a[k] - mx);
        ll += mx + std::log(den);
        for (int k = 0; k < K; ++k)
            r[static_cast<std::size_t>(i) * K + k] = std::exp(a[k] - mx) / den;
    }
    return static_cast<double>(ll);
}

struct Point {
    double t, mu, f, dx;
    Point(double t_ = 0, double m_ = 0, double f_ = INF, double dx_ = 0)
        : t(t_), mu(m_), f(f_), dx(dx_) {}
};

// The responsibility weights are normalized within each component. This
// leaves its optimum unchanged and protects very small component masses.
struct Profile {
    const Data& d;
    std::vector<int> idx;
    std::vector<double> a, loga;
    double xmax;
    Profile(const Data& data, const std::vector<double>& r, int K, int k)
        : d(data), xmax(0) {
        long double mass = 0;
        for (int i = 0; i < d.n; ++i) mass += r[static_cast<std::size_t>(i) * K + k];
        if (mass <= 0) throw std::runtime_error("Cannot profile an empty component.");
        for (int i = 0; i < d.n; ++i) {
            double ri = r[static_cast<std::size_t>(i) * K + k];
            if (ri > 0) {
                idx.push_back(i);
                a.push_back(static_cast<double>(ri / mass));
                loga.push_back(std::log(ri) - static_cast<double>(std::log(mass)));
            }
        }
        xmax = x_of_t(d.upper);
    }
    double x_of_t(double t) const {
        if (t <= 0) return 0;
        double u = t / d.vmin;
        if (std::isfinite(u)) return std::log1p(u);
        return std::log(t) - std::log(d.vmin) + std::log1p(d.vmin / t);
    }
    double t_of_x(double x) const {
        if (x <= 0) return 0;
        if (x >= xmax) return d.upper;
        double t = (x < 35) ? d.vmin * std::expm1(x) :
            std::exp(std::log(d.vmin) + x) - d.vmin;
        return clamp(t, 0, d.upper);
    }
    Point at_t(double t) const {
        // Conditional mean computed with rescaled inverse-variance weights.
        double maxlogw = -INF;
        for (std::size_t j = 0; j < idx.size(); ++j)
            maxlogw = std::max(maxlogw, loga[j] - std::log(d.v[idx[j]] + t));
        long double sw = 0, swy = 0;
        for (std::size_t j = 0; j < idx.size(); ++j) {
            double w = std::exp(loga[j] - std::log(d.v[idx[j]] + t) - maxlogw);
            sw += w;
            swy += static_cast<long double>(w) * d.y[idx[j]];
        }
        double mu = static_cast<double>(swy / sw);
        long double f = 0, dx = 0;
        for (std::size_t j = 0; j < idx.size(); ++j) {
            double s = d.v[idx[j]] + t;
            long double z = static_cast<long double>(d.y[idx[j]]) - mu;
            long double z2s = z * z / s;
            f += static_cast<long double>(a[j]) * (std::log(s) + z2s);
            // Envelope theorem; derivative with respect to x = log1p(t/vmin).
            dx += static_cast<long double>(a[j]) * ((d.vmin + t) / s) * (1 - z2s);
        }
        return Point(t, mu, static_cast<double>(0.5L * f), static_cast<double>(0.5L * dx));
    }
    Point at_x(double x) const { return at_t(t_of_x(x)); }
};

inline void consider(Point& best, const Point& p) {
    if (std::isfinite(p.f) &&
        (p.f < best.f || (p.f == best.f && p.t < best.t))) best = p;
}

// Bounded Brent search within one candidate basin. Global robustness comes
// from scanning and refining multiple basins, not assuming unimodality.
Point brent(const Profile& f, double lo, double hi, double x, double tol) {
    const double golden = 0.3819660112501051518;
    double w = x, v = x, e = 0, step = 0;
    Point px = f.at_x(x), best = px;
    double fx = px.f, fw = fx, fv = fx;
    consider(best, f.at_x(lo)); consider(best, f.at_x(hi));
    for (int it = 0; it < 160; ++it) {
        double mid = 0.5 * (lo + hi);
        double eps = tol * (1 + std::fabs(x));
        double eps2 = 2 * eps;
        if (std::fabs(x - mid) <= eps2 - 0.5 * (hi - lo)) break;
        bool parabolic = false;
        if (std::fabs(e) > eps && std::isfinite(fx) && std::isfinite(fw) && std::isfinite(fv)) {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2 * (q - r);
            if (q > 0) p = -p;
            q = std::fabs(q);
            double olde = e; e = step;
            if (q > 0 && std::fabs(p) < std::fabs(0.5 * q * olde) &&
                p > q * (lo - x) && p < q * (hi - x)) {
                step = p / q;
                double u = x + step;
                if (u - lo < eps2 || hi - u < eps2) step = (mid >= x) ? eps : -eps;
                parabolic = true;
            }
        }
        if (!parabolic) { e = (x < mid) ? hi - x : lo - x; step = golden * e; }
        double u = x + ((std::fabs(step) >= eps) ? step : ((step >= 0) ? eps : -eps));
        u = clamp(u, lo, hi);
        Point pu = f.at_x(u); consider(best, pu);
        double fu = pu.f;
        if (fu <= fx) {
            if (u >= x) lo = x; else hi = x;
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        } else {
            if (u < x) lo = u; else hi = u;
            if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    return best;
}

Point optimize_profile(const Profile& f, double current, int grid_size, double tol) {
    Point best = f.at_t(current); // Never discard the previous variance.
    consider(best, f.at_t(0));
    if (f.d.upper == 0) return best;
    consider(best, f.at_t(f.d.upper));
    std::vector<double> xs;
    xs.reserve(2 * grid_size + f.d.vknots.size() + 1);
    xs.push_back(f.x_of_t(current));
    for (int j = 0; j < grid_size; ++j) {
        double u = static_cast<double>(j) / (grid_size - 1);
        xs.push_back(u * f.xmax);                    // Logarithmic variance grid.
        xs.push_back(f.x_of_t(u * f.d.upper));        // Linear variance grid.
    }
    for (double t : f.d.vknots) xs.push_back(f.x_of_t(t));
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    std::vector<Point> val;
    val.reserve(xs.size());
    for (double x : xs) { val.push_back(f.at_x(x)); consider(best, val.back()); }
    // Refine every sampled local minimum.
    for (std::size_t j = 1; j + 1 < xs.size(); ++j)
        if (val[j].f <= val[j-1].f && val[j].f <= val[j+1].f &&
            (val[j].f < val[j-1].f || val[j].f < val[j+1].f))
            consider(best, brent(f, xs[j-1], xs[j+1], xs[j], tol));
    // Also refine every bracketed negative-to-positive profile-score crossing.
    for (std::size_t j = 0; j + 1 < xs.size(); ++j) {
        if (!(val[j].dx < 0 && val[j+1].dx > 0)) continue;
        double a = xs[j], b = xs[j+1];
        for (int it = 0; it < 100; ++it) {
            double m = 0.5 * (a + b);
            Point pm = f.at_x(m); consider(best, pm);
            if (b - a <= tol * (1 + std::fabs(m)) || pm.dx == 0) break;
            if (pm.dx < 0) a = m; else b = m;
        }
    }
    return best;
}

Params initialize(const Data& d, int K, int start) {
    Params p(K);
    std::vector<double> sy = d.y;
    std::sort(sy.begin(), sy.end());
    long double av = 0, avv = 0, var = 0;
    for (int i = 0; i < d.n; ++i) { av += d.y[i]; avv += d.v[i]; }
    av /= d.n; avv /= d.n;
    for (double y : d.y) var += (y - av) * (y - av);
    var /= d.n;
    double moment = std::max(0.0, static_cast<double>(var - avv));
    double sump = 0;
    for (int k = 0; k < K; ++k) {
        if (start == 0) {
            int pos = std::min(d.n - 1, static_cast<int>((k + 0.5) * d.n / K));
            p.mu[k] = sy[pos];
            p.tau[k] = std::max(moment, d.upper / (16.0 * K * K));
            p.pi[k] = 1.0;
        } else {
            // R::runif draws from R's current RNG stream.
            int index = std::min(d.n - 1,
                static_cast<int>(R::runif(0.0, static_cast<double>(d.n))));
            p.mu[k] = d.y[index];
            if (d.upper == 0) p.tau[k] = 0;
            else {
                double big = d.upper * 0.5;
                if (big == 0) big = d.upper;
                double small = std::max(d.upper * 1e-8, std::min(d.vmin, big * 0.5));
                if (small == 0) small = big;
                p.tau[k] = std::exp(std::log(small) +
                    R::runif(0.0, 1.0) * (std::log(big) - std::log(small)));
            }
            p.pi[k] = 0.25 + R::runif(0.0, 1.0);
        }
        p.tau[k] = clamp(p.tau[k], 0, d.upper);
        sump += p.pi[k];
    }
    for (double& q : p.pi) q /= sump;
    return p;
}


// Equality constraints are represented by a reduced parameter vector.
struct Constraint {
    bool equal_mu, equal_tau;
    int K, nm, nt, q;
    Constraint(int k, bool em, bool et): equal_mu(em), equal_tau(et), K(k),
        nm(em?1:k), nt(et?1:k), q(nm+nt+k-1) {}
    void enforce(Params& p) const {
        if (equal_mu) {
            double avg=0; for (int k=0;k<K;++k) avg+=p.pi[k]*p.mu[k];
            std::fill(p.mu.begin(),p.mu.end(),avg);
        }
        if (equal_tau) {
            double avg=0; for (int k=0;k<K;++k) avg+=p.pi[k]*p.tau[k];
            std::fill(p.tau.begin(),p.tau.end(),avg);
        }
    }
    std::vector<double> pack(const Params& p) const {
        std::vector<double> th(q);
        for (int j=0;j<nm;++j) th[j]=p.mu[j];
        for (int j=0;j<nt;++j) th[nm+j]=p.tau[j];
        for (int j=0;j<K-1;++j)
            th[nm+nt+j]=std::log(std::max(1e-12,p.pi[j]))-std::log(std::max(1e-12,p.pi[K-1]));
        return th;
    }
    Params unpack(const std::vector<double>& th, std::vector<double>& lp) const {
        if ((int)th.size()!=q) throw std::invalid_argument("Incorrect constrained parameter length");
        Params p(K);
        double mx=0;
        for(int j=0;j<K-1;++j) mx=std::max(mx,th[nm+nt+j]);
        double den=std::exp(-mx);
        for(int j=0;j<K-1;++j) den+=std::exp(th[nm+nt+j]-mx);
        double logden=mx+std::log(den);
        lp.resize(K);
        for(int k=0;k<K;++k){
            p.mu[k]=th[equal_mu?0:k];
            p.tau[k]=th[nm+(equal_tau?0:k)];
            lp[k]=(k<K-1?th[nm+nt+k]:0)-logden;
            p.pi[k]=std::exp(lp[k]);
        }
        return p;
    }
};

// The same reduced objective/gradient is used for refinement and SEs.
struct ConstrainedJoint {
    Data d;
    Constraint c;
    std::vector<double> last, gradient;
    double value;
    bool cached;
    ConstrainedJoint(const Data& data, const Constraint& cc):d(data),c(cc),value(INF),cached(false){}
    Params unpack(const std::vector<double>& th,std::vector<double>& lp) const {return c.unpack(th,lp);}
    void evaluate(const std::vector<double>& th){
        if(cached && th==last)return;
        for(double z:th)if(!std::isfinite(z))throw std::runtime_error("Nonfinite constrained parameter");
        std::vector<double> lp,r;
        Params p=c.unpack(th,lp);
        for(double t:p.tau)if(t<0)throw std::runtime_error("Negative variance");
        double ll=estep(d,p,r,&lp);
        if(!std::isfinite(ll))throw std::runtime_error("Nonfinite likelihood");
        std::vector<long double> g(c.q,0);
        for(int i=0;i<d.n;++i)for(int k=0;k<c.K;++k){
            double rik=r[(std::size_t)i*c.K+k];
            if(rik==0)continue;
            long double v=d.v[i]+p.tau[k], z=(long double)p.mu[k]-d.y[i];
            g[c.equal_mu?0:k]+=rik*z/v;
            g[c.nm+(c.equal_tau?0:k)]+=0.5L*rik*(1-z*z/v)/v;
            if(k<c.K-1)g[c.nm+c.nt+k]-=rik;
        }
        gradient.resize(c.q);
        for(int j=0;j<c.q;++j)gradient[j]=(double)(g[j]/d.n);
        for(int j=0;j<c.K-1;++j)gradient[c.nm+c.nt+j]+=p.pi[j];
        value=-ll/d.n;last=th;cached=true;
    }
};

// Robust bounded one-dimensional search: dense linear/log grids and Brent
// refinement around every sampled local minimum.
template<class F> double bounded_minimum(F fun,double current,double upper,int grid_size,double vmin,double tol){
    if(upper<=0)return 0;
    std::vector<double> xs;
    xs.push_back(clamp(current,0,upper));
    for(int j=0;j<grid_size;++j){
        double u=(double)j/(grid_size-1);
        xs.push_back(u*upper);
        xs.push_back(vmin*std::expm1(u*std::log1p(upper/vmin)));
    }
    std::sort(xs.begin(),xs.end());
    xs.erase(std::unique(xs.begin(),xs.end()),xs.end());
    std::vector<double> fs(xs.size());
    int best=0;
    for(std::size_t j=0;j<xs.size();++j){fs[j]=fun(xs[j]);if(fs[j]<fs[best])best=(int)j;}
    double bt=xs[best],bf=fs[best];
    for(std::size_t j=1;j+1<xs.size();++j){
        if(fs[j]>fs[j-1]||fs[j]>fs[j+1])continue;
        double a=xs[j-1],b=xs[j+1];
        const double phi=0.3819660112501051;
        double x1=a+phi*(b-a),x2=b-phi*(b-a),f1=fun(x1),f2=fun(x2);
        for(int it=0;it<130 && b-a>tol*(1+std::fabs(a)+std::fabs(b));++it){
            if(f1<f2){b=x2;x2=x1;f2=f1;x1=a+phi*(b-a);f1=fun(x1);}
            else {a=x1;x1=x2;f1=f2;x2=b-phi*(b-a);f2=fun(x2);}
        }
        if(f1<bf){bf=f1;bt=x1;}
        if(f2<bf){bf=f2;bt=x2;}
    }
    return bt;
}

// Constrained generalized ECM: exact conditional mean updates, followed by
// variance coordinate updates. Every outer step is likelihood-safeguarded.
struct Run {
    Params p;
    double ll;
    int iterations;
    bool converged;
    std::string status;
    std::vector<double> trace;
    Run() : ll(-INF), iterations(0), converged(false), status("not_run") {}
};

Run ecm(const Data& d, Params p, int maxit, double tol,
        int grid_size, double profile_tol, bool keep_trace) {
    Run out;
    const int K = p.size();
    std::vector<double> r;
    double ll = estep(d, p, r);
    out.p = p; out.ll = ll;
    if (!std::isfinite(ll)) { out.status = "nonfinite_start"; return out; }
    if (keep_trace) out.trace.push_back(ll);
    out.status = "maxit";
    for (int it = 0; it < maxit; ++it) {
        if (it % 5 == 0) interrupt();
        Params next = p;
        // All component updates use the SAME responsibilities from the E-step.
        for (int k = 0; k < K; ++k) {
            long double mass = 0;
            for (int i = 0; i < d.n; ++i) mass += r[static_cast<std::size_t>(i) * K + k];
            next.pi[k] = static_cast<double>(mass / d.n);
            if (mass > 0) {
                Profile f(d, r, K, k);
                Point opt = optimize_profile(f, p.tau[k], grid_size, profile_tol);
                next.mu[k] = opt.mu;
                next.tau[k] = opt.t;
            } // Exactly empty components retain their unidentified mu and tau.
        }
        double sump = std::accumulate(next.pi.begin(), next.pi.end(), 0.0);
        for (double& q : next.pi) q /= sump; // Floating-point normalization only.
        std::vector<double> nr;
        double nll = estep(d, next, nr);
        if (!std::isfinite(nll)) { out.status = "nonfinite_update_rejected"; break; }
        if (nll < ll) {
            double roundoff = 128 * std::numeric_limits<double>::epsilon() * (1 + std::fabs(ll));
            out.converged = (ll - nll <= roundoff);
            out.status = out.converged ? "roundoff_limit" : "likelihood_decrease_rejected";
            break; // Retain the old parameters: stored likelihood never decreases.
        }
        double gain_per_study = (nll - ll) / d.n;
        p = std::move(next); r = std::move(nr); ll = nll;
        out.iterations = it + 1;
        if (keep_trace) out.trace.push_back(ll);
        if (gain_per_study <= tol) {
            out.converged = true; out.status = "likelihood_tolerance"; break;
        }
    }
    out.p = p; out.ll = ll;
    return out;
}


Run ecm_constrained(const Data& d,Params p,const Constraint& c,int maxit,double tol,
                    int grid_size,double profile_tol,bool keep_trace){
    Run out;
    c.enforce(p);
    std::vector<double> r;
    double ll=estep(d,p,r);
    out.p=p;out.ll=ll;
    if(!std::isfinite(ll)){out.status="nonfinite_start";return out;}
    if(keep_trace)out.trace.push_back(ll);
    out.status="maxit";
    for(int it=0;it<maxit;++it){
        if(it%5==0)interrupt();
        Params next=p;
        for(int k=0;k<c.K;++k){
            long double mass=0;
            for(int i=0;i<d.n;++i)mass+=r[(std::size_t)i*c.K+k];
            next.pi[k]=(double)(mass/d.n);
        }
        double psum=std::accumulate(next.pi.begin(),next.pi.end(),0.0);
        for(double& w:next.pi)w/=psum;
        // Inner coordinate ascent uses fixed responsibilities.
        for(int inner=0;inner<12;++inner){
            Params old=next;
            // Conditional means (shared or class-specific).
            if(c.equal_mu){
                long double sw=0,swy=0;
                for(int i=0;i<d.n;++i)for(int k=0;k<c.K;++k){
                    long double w=r[(std::size_t)i*c.K+k]/(d.v[i]+next.tau[k]);
                    sw+=w;swy+=w*d.y[i];
                }
                if(sw>0)std::fill(next.mu.begin(),next.mu.end(),(double)(swy/sw));
            }else{
                for(int k=0;k<c.K;++k){
                    long double sw=0,swy=0;
                    for(int i=0;i<d.n;++i){
                        long double w=r[(std::size_t)i*c.K+k]/(d.v[i]+next.tau[k]);
                        sw+=w;swy+=w*d.y[i];
                    }
                    if(sw>0)next.mu[k]=(double)(swy/sw);
                }
            }
            if(c.equal_tau){
                auto fun=[&](double t){
                    long double f=0;
                    for(int i=0;i<d.n;++i)for(int k=0;k<c.K;++k){
                        long double w=r[(std::size_t)i*c.K+k],v=d.v[i]+t,z=d.y[i]-next.mu[k];
                        f+=w*(std::log((double)v)+z*z/v);
                    }
                    return (double)(f/(2*d.n));
                };
                double t=bounded_minimum(fun,next.tau[0],d.upper,grid_size,d.vmin,profile_tol);
                std::fill(next.tau.begin(),next.tau.end(),t);
            }else{
                for(int k=0;k<c.K;++k){
                    auto fun=[&](double t){
                        long double f=0;
                        for(int i=0;i<d.n;++i){
                            long double w=r[(std::size_t)i*c.K+k],v=d.v[i]+t,z=d.y[i]-next.mu[k];
                            f+=w*(std::log((double)v)+z*z/v);
                        }
                        return (double)(f/d.n);
                    };
                    next.tau[k]=bounded_minimum(fun,next.tau[k],d.upper,grid_size,d.vmin,profile_tol);
                }
            }
            double change=0;
            for(int k=0;k<c.K;++k)change=std::max(change,std::max(std::fabs(next.mu[k]-old.mu[k]),std::fabs(next.tau[k]-old.tau[k])));
            if(change<1e-10)break;
        }
        std::vector<double> nr;
        double nll=estep(d,next,nr);
        if(!std::isfinite(nll)||nll<ll){
            double roundoff=128*std::numeric_limits<double>::epsilon()*(1+std::fabs(ll));
            out.converged=std::isfinite(nll)&&ll-nll<=roundoff;
            out.status=out.converged?"roundoff_limit":"likelihood_decrease_rejected";
            break;
        }
        double gain=(nll-ll)/d.n;
        p=std::move(next);r=std::move(nr);ll=nll;
        out.iterations=it+1;
        if(keep_trace)out.trace.push_back(ll);
        if(gain<=tol){out.converged=true;out.status="likelihood_tolerance";break;}
    }
    out.p=p;out.ll=ll;
    return out;
}

// Parameter order for joint optimization:
// (mu[1:K], tau2[1:K], eta[1:(K-1)]), eta[k] = log(pi[k]/pi[K]).
// Objective is minus the observed log-likelihood PER STUDY.
struct Joint {
    Data d;
    int K;
    std::vector<double> last, gradient;
    double value;
    bool cached;
    Joint(const Data& data, int k) : d(data), K(k), value(INF), cached(false) {}
    Params unpack(const std::vector<double>& theta, std::vector<double>& logpi) const {
        if (theta.size() != static_cast<std::size_t>(3 * K - 1))
            throw std::invalid_argument("Incorrect joint parameter length.");
        Params p(K);
        double mx = 0;
        for (int k = 0; k < K - 1; ++k) mx = std::max(mx, theta[2*K+k]);
        double den = std::exp(-mx);
        for (int k = 0; k < K - 1; ++k) den += std::exp(theta[2*K+k] - mx);
        double logden = mx + std::log(den);
        logpi.resize(K);
        for (int k = 0; k < K; ++k) {
            p.mu[k] = theta[k]; p.tau[k] = theta[K+k];
            logpi[k] = ((k < K-1) ? theta[2*K+k] : 0.0) - logden;
            p.pi[k] = std::exp(logpi[k]);
        }
        return p;
    }
    void evaluate(const std::vector<double>& theta) {
        if (cached && theta == last) return;
        for (double x : theta)
            if (!std::isfinite(x)) throw std::runtime_error("Nonfinite joint parameter.");
        std::vector<double> lp, r;
        Params p = unpack(theta, lp);
        for (int k = 0; k < K; ++k)
            if (p.tau[k] < 0) throw std::runtime_error("Negative joint variance.");
        double ll = estep(d, p, r, &lp);
        if (!std::isfinite(ll)) throw std::runtime_error("Nonfinite joint likelihood.");
        std::vector<long double> g(theta.size(), 0);
        for (int i = 0; i < d.n; ++i) {
            for (int k = 0; k < K; ++k) {
                double rik = r[static_cast<std::size_t>(i) * K + k];
                if (rik == 0) continue;
                long double s = d.v[i] + p.tau[k];
                long double z = static_cast<long double>(p.mu[k]) - d.y[i];
                g[k] += rik * z / s;
                g[K+k] += 0.5L * rik * (1 - z*z/s) / s;
                if (k < K-1) g[2*K+k] -= rik;
            }
        }
        gradient.resize(theta.size());
        for (std::size_t j = 0; j < g.size(); ++j)
            gradient[j] = static_cast<double>(g[j] / d.n);
        for (int k = 0; k < K-1; ++k) gradient[2*K+k] += p.pi[k];
        for (double x : gradient)
            if (!std::isfinite(x)) throw std::runtime_error("Nonfinite joint gradient.");
        value = -ll / d.n; last = theta; cached = true;
    }
};

Params sorted(const Params& p) {
    int K = p.size();
    std::vector<int> ix(K); std::iota(ix.begin(), ix.end(), 0);
    std::stable_sort(ix.begin(), ix.end(), [&](int a, int b) { return p.mu[a] < p.mu[b]; });
    Params ans(K);
    for (int k = 0; k < K; ++k) {
        ans.mu[k] = p.mu[ix[k]]; ans.tau[k] = p.tau[ix[k]]; ans.pi[k] = p.pi[ix[k]];
    }
    return ans;
}

} // namespace metamixture_core

namespace metamixture_rcpp {

double joint_fn(Rcpp::NumericVector par, SEXP context) {
    Rcpp::XPtr<metamixture_core::Joint> obj(context);
    obj->evaluate(Rcpp::as<std::vector<double> >(par));
    return obj->value;
}
Rcpp::NumericVector joint_gr(Rcpp::NumericVector par, SEXP context) {
    Rcpp::XPtr<metamixture_core::Joint> obj(context);
    obj->evaluate(Rcpp::as<std::vector<double> >(par));
    return Rcpp::wrap(obj->gradient);
}
Rcpp::List params_list(const metamixture_core::Data& d, const metamixture_core::Params& p) {
    int K = p.size();
    Rcpp::NumericVector mu(K), tau(K), pi(K);
    for (int k = 0; k < K; ++k) {
        mu[k] = d.center + d.scale * p.mu[k];
        tau[k] = (p.tau[k] * d.scale) * d.scale;
        pi[k] = p.pi[k];
        if (!std::isfinite(mu[k]) || !std::isfinite(tau[k]))
            Rcpp::stop("An estimate overflows on the original measurement scale.");
    }
    return Rcpp::List::create(Rcpp::_["mu"] = mu, Rcpp::_["tau2"] = tau, Rcpp::_["pi"] = pi);
}

double constrained_fn(Rcpp::NumericVector par,SEXP context){
    Rcpp::XPtr<metamixture_core::ConstrainedJoint> obj(context);
    obj->evaluate(Rcpp::as<std::vector<double> >(par));return obj->value;
}
Rcpp::NumericVector constrained_gr(Rcpp::NumericVector par,SEXP context){
    Rcpp::XPtr<metamixture_core::ConstrainedJoint> obj(context);
    obj->evaluate(Rcpp::as<std::vector<double> >(par));return Rcpp::wrap(obj->gradient);
}

metamixture_core::Params read_init(Rcpp::List init, const metamixture_core::Data& d, int K) {
    if (!init.containsElementNamed("mu") || !init.containsElementNamed("tau2") ||
        !init.containsElementNamed("pi")) Rcpp::stop("init must contain mu, tau2, and pi.");
    std::vector<double> mu = Rcpp::as<std::vector<double> >(init["mu"]);
    std::vector<double> tau = Rcpp::as<std::vector<double> >(init["tau2"]);
    std::vector<double> pi = Rcpp::as<std::vector<double> >(init["pi"]);
    if (mu.size() != static_cast<std::size_t>(K) || tau.size() != mu.size() || pi.size() != mu.size())
        Rcpp::stop("Each vector in init must have length K.");
    metamixture_core::Params p(K);
    for (int k = 0; k < K; ++k) {
        if (!std::isfinite(mu[k]) || !std::isfinite(tau[k]) || tau[k] < 0 ||
            !std::isfinite(pi[k]) || pi[k] < 0)
            Rcpp::stop("init requires finite mu, nonnegative tau2, and nonnegative pi.");
        double m = (mu[k] - d.center) / d.scale;
        double t = (tau[k] / d.scale) / d.scale;
        if (!std::isfinite(m) || !std::isfinite(t)) Rcpp::stop("init cannot be scaled safely.");
        // Restrict starts to a box that contains a global maximum.
        p.mu[k] = metamixture_core::clamp(m, d.ymin, d.ymax);
        p.tau[k] = metamixture_core::clamp(t, 0, d.upper);
    }
    double largest = *std::max_element(pi.begin(), pi.end());
    if (largest <= 0) Rcpp::stop("At least one initial mixing proportion must be positive.");
    long double total = 0;
    for (double w : pi) total += w / largest;
    for (int k = 0; k < K; ++k) p.pi[k] = static_cast<double>((pi[k] / largest) / total);
    return p;
}

// Observed-information standard errors. The core already supplies the exact
// analytic gradient of -logLik/n on internally standardized data.
// Central differences of that gradient yield the observed Hessian.
Rcpp::List standard_errors(const metamixture_core::Data& d,
                           const metamixture_core::Params& p,
                           const std::vector<double>& yi,
                           double rel_step, double boundary_tol,
                           bool converged, const metamixture_core::Constraint& c) {
    (void)converged;
    using namespace Rcpp;
    const int K = p.size(), q = c.q, m = 3*K;
    if (!std::isfinite(rel_step) || rel_step <= 0 || rel_step >= .1)
        stop("se_rel_step must be in (0, 0.1).");
    if (!std::isfinite(boundary_tol) || boundary_tol < 0)
        stop("se_boundary_tol must be finite and nonnegative.");
    double mean = std::accumulate(yi.begin(), yi.end(), 0.0)/yi.size();
    long double ss = 0;
    for (double z : yi) ss += (static_cast<long double>(z)-mean)*(z-mean);
    double variance = static_cast<double>(ss/(yi.size()-1));
    if (!std::isfinite(variance)) throw std::runtime_error("Cannot compute sample variance for SE boundary check.");
    for (int k=0;k<K;++k) {
        const double tau_orig = (p.tau[k]*d.scale)*d.scale;
        if (tau_orig <= boundary_tol*std::max(1.0, variance) ||
            p.pi[k] <= boundary_tol || (K>1 && p.pi[k]>=1-boundary_tol))
            throw std::runtime_error("Boundary fit: inverse-information standard errors are unreliable.");
    }
    std::vector<double> theta=c.pack(p);
    metamixture_core::ConstrainedJoint joint(d,c);
    NumericMatrix H(q,q);
    for (int j=0;j<q;++j) {
        double h=rel_step*std::max(1.0,std::fabs(theta[j]));
        if (j>=c.nm && j<c.nm+c.nt) h=std::min(h, theta[j]*.25);
        if (!(h>0) || !std::isfinite(h))
            throw std::runtime_error("Unable to construct interior finite-difference step.");
        std::vector<double> plus=theta, minus=theta;
        plus[j]+=h; minus[j]-=h;
        joint.evaluate(plus);
        std::vector<double> gp=joint.gradient;
        joint.evaluate(minus);
        std::vector<double> gm=joint.gradient;
        for (int i=0;i<q;++i)
            H(i,j) = -d.n*(gp[i]-gm[i])/(2*h);
    }
    NumericMatrix info(q,q);
    for (int j=0;j<q;++j)
        for (int i=0;i<q;++i)
            info(i,j)=-.5*(H(i,j)+H(j,i));
    Environment base=Environment::base_env();
    Function eigen=base["eigen"], solve=base["solve"];
    List eig=eigen(info, _["symmetric"]=true, _["only.values"]=true);
    NumericVector ev=eig["values"];
    double lo=*std::min_element(ev.begin(),ev.end());
    double hi=*std::max_element(ev.begin(),ev.end());
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo<=0 || lo<=hi*1e-8)
        throw std::runtime_error("Observed information is singular or nearly singular.");
    NumericMatrix V=solve(info);
    NumericMatrix J(m,q);
    for (int k=0;k<K;++k) {
        J(k,c.equal_mu?0:k)=d.scale;
        J(K+k,c.nm+(c.equal_tau?0:k))=d.scale*d.scale;
        for (int j=0;j<K-1;++j)
            J(2*K+k,c.nm+c.nt+j)=p.pi[k]*((k==j?1.0:0.0)-p.pi[j]);
    }
    NumericMatrix Vout(m,m);
    // Compute J V J^T without external linear-algebra dependencies.
    std::vector<double> JV(static_cast<std::size_t>(m)*q,0.0);
    for (int a=0;a<m;++a) for (int j=0;j<q;++j)
        for (int b=0;b<q;++b) JV[a*q+j]+=J(a,b)*V(b,j);
    for (int a=0;a<m;++a) for (int b=0;b<m;++b)
        for (int j=0;j<q;++j) Vout(a,b)+=JV[a*q+j]*J(b,j);
    CharacterVector nm(m);
    NumericVector estimates(m), ses(m), se_mu(K), se_tau2(K), se_pi(K);
    for (int k=0;k<K;++k) {
        nm[k]="mu"+std::to_string(k+1);
        nm[K+k]="tau2"+std::to_string(k+1);
        nm[2*K+k]="pi"+std::to_string(k+1);
        estimates[k]=d.center+d.scale*p.mu[k];
        estimates[K+k]=(p.tau[k]*d.scale)*d.scale;
        estimates[2*K+k]=p.pi[k];
    }
    for (int a=0;a<m;++a) {
        if (!std::isfinite(Vout(a,a))) throw std::runtime_error("Nonfinite SE covariance diagonal.");
        ses[a]=std::sqrt(std::max(0.0,Vout(a,a)));
    }
    for (int k=0;k<K;++k) {
        se_mu[k]=ses[k]; se_tau2[k]=ses[K+k]; se_pi[k]=ses[2*K+k];
    }
    ses.attr("names")=nm; estimates.attr("names")=nm;
    Vout.attr("dimnames")=List::create(nm,nm);
    return List::create(_["se_mu"]=se_mu, _["se_tau2"]=se_tau2,
        _["se_pi"]=se_pi, _["estimates"]=estimates, _["se"]=ses,
        _["vcov"]=Vout, _["information"]=info,
        _["min_information_eigenvalue"]=lo,
        _["method"]="Inverse observed information; central differences of analytic score");
}

} // namespace metamixture_rcpp

//' Fit a random-effects mixture model
//'
//' Maximum-likelihood estimation of a finite mixture of normal random-effects
//' meta-analysis models with known sampling variances. Uses safeguarded
//' profiled component updates and optional joint L-BFGS-B refinement.
//'
//' @param yi Numeric vector of observed effect sizes; all values must be finite.
//' @param vi Numeric vector of known sampling variances, of the same length as
//'   yi; all values must be finite and strictly positive.
//' @param K Number of components, between 1 and length(yi).
//' @param nstart Number of starting solutions. The first is deterministic
//'   unless replaced by init; later starts use R's random number generator.
//' @param maxit Maximum number of ECM iterations per start.
//' @param tol ECM stopping tolerance on log-likelihood gain per observation.
//' @param grid_size Number of points in each component-profile search grid.
//' @param profile_tol Local search tolerance in the transformed coordinate
//'   log1p(tau2 / min(vi)), using internally scaled data.
//' @param joint Whether to refine the best ECM fit jointly with L-BFGS-B.
//' @param joint_maxit Maximum number of joint optimizer iterations.
//' @param joint_factr L-BFGS-B relative objective tolerance control.
//' @param joint_pgtol L-BFGS-B projected-gradient tolerance control.
//' @param init Optional list with length-K vectors mu, tau2, and pi. Replaces
//'   the first starting solution rather than adding an extra start.
//' @param keep_trace Whether to return the best ECM start's likelihood trace.
//' @param compute_se Calculate observed-information standard errors (default TRUE).
//' @param se_rel_step Relative step for numerical differentiation of the analytic score.
//' @param se_boundary_tol Boundary tolerance; SEs are NA when regular Wald inference is inappropriate.
//' @param mu Whether component means are `varying` or `equal`.
//' @param tau2 Whether component heterogeneity variances are `varying` or `equal`.
//'
//' @return A list with component estimates mu, tau2, and pi; posterior
//'   probabilities and classifications; effective_n; logLik, AIC, and BIC;
//'   n, K, and npar; converged; ecm and joint diagnostics; starts; and notes.
//'   Equality-constrained estimates are repeated for all K components.
//'   When requested, also returns se_mu, se_tau2, se_pi, se, vcov, information,
//'   se_details, and se_status. On a boundary or singular fit, SEs are NA.
//'   Components are sorted by increasing estimated mean. The ECM and joint
//'   convergence diagnostics should be inspected separately.
//'
//' @details Random initializations respect R's RNG state. Use set.seed() in R
//'   before calling this function for reproducibility. There is no seed
//'   argument and the function does not seed a private generator. The Rcpp
//'   attributes wrapper manages RNG state using an RNGScope. Multiple starts
//'   do not guarantee a global maximum. Estimation is ML, not REML.
//'
//' @useDynLib metamixture, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
Rcpp::List remix_cpp(
        Rcpp::NumericVector yi, Rcpp::NumericVector vi,
        int K = 2, int nstart = 20, int maxit = 1000,
        double tol = 1e-9, int grid_size = 65, double profile_tol = 1e-10,
        bool joint = true, int joint_maxit = 1000,
        double joint_factr = 1e7, double joint_pgtol = 1e-8,
        Rcpp::Nullable<Rcpp::List> init = R_NilValue,
        bool keep_trace = false, bool compute_se = true,
        double se_rel_step = 1e-4, double se_boundary_tol = 1e-10,
        std::string mu = "varying", std::string tau2 = "varying") {
    using namespace Rcpp;
    if (K == NA_INTEGER || nstart == NA_INTEGER || maxit == NA_INTEGER ||
        grid_size == NA_INTEGER || joint_maxit == NA_INTEGER ||
        K < 1 || K > yi.size() || nstart < 1 || maxit < 1 || grid_size < 9 ||
        joint_maxit < 1 || !std::isfinite(tol) || tol <= 0 ||
        !std::isfinite(profile_tol) || profile_tol <= 0 || profile_tol > 0.1 ||
        !std::isfinite(joint_factr) || joint_factr <= 0 ||
        !std::isfinite(joint_pgtol) || joint_pgtol < 0)
        stop("Invalid K or optimization controls.");
    if (mu != "varying" && mu != "equal") stop("mu must be varying or equal");
    if (tau2 != "varying" && tau2 != "equal") stop("tau2 must be varying or equal");
    if (K > 1 && mu == "equal" && tau2 == "equal")
        stop("When both mu and tau2 are equal, the K classes are identical and mixing proportions are unidentified. Fit K=1 instead.");
    metamixture_core::Constraint constraints(K, mu == "equal", tau2 == "equal");
    metamixture_core::Data d(as<std::vector<double> >(yi), as<std::vector<double> >(vi));
    metamixture_core::Params user_init;
    bool have_init = init.isNotNull();
    if (have_init) { user_init = metamixture_rcpp::read_init(List(init.get()), d, K); constraints.enforce(user_init); }
    metamixture_core::Run best;
    int best_start = -1;
    NumericVector start_ll(nstart);
    IntegerVector start_it(nstart), active(nstart), start_number(nstart);
    LogicalVector start_conv(nstart);
    CharacterVector start_status(nstart);
    for (int s = 0; s < nstart; ++s) {
        checkUserInterrupt();
        metamixture_core::Params p = (s == 0 && have_init) ? user_init : metamixture_core::initialize(d, K, s);
        constraints.enforce(p);
        metamixture_core::Run run = (!constraints.equal_mu && !constraints.equal_tau) ?
            metamixture_core::ecm(d, p, maxit, tol, grid_size, profile_tol, keep_trace) :
            metamixture_core::ecm_constrained(d, p, constraints, maxit, tol, grid_size, profile_tol, keep_trace);
        start_number[s] = s + 1; start_ll[s] = d.original_ll(run.ll);
        start_it[s] = run.iterations; start_conv[s] = run.converged;
        start_status[s] = run.status;
        int a = 0; for (double w : run.p.pi) if (w > 0) ++a;
        active[s] = a;
        if (run.ll > best.ll ||
            (run.ll == best.ll && run.converged && !best.converged)) {
            best = std::move(run); best_start = s + 1;
        }
    }
    if (best_start < 0) stop("No finite ECM solution was found. Check input scales and initialization.");
    metamixture_core::Params final_p = best.p;
    bool joint_attempted = false, joint_accepted = false;
    bool joint_start_adjusted = false;
    int joint_code = NA_INTEGER;
    double joint_ll = NA_REAL, joint_start_ll = NA_REAL, pg_norm = NA_REAL;
    std::string joint_message = joint ? "not_run" : "disabled";
    IntegerVector joint_counts = IntegerVector::create(_["function"] = NA_INTEGER, _["gradient"] = NA_INTEGER);

    if (joint && d.upper == 0) {
        joint_message = "Skipped: all yi are equal; means and variances are fixed at their optimum.";
    } else if (joint) {
        joint_attempted = true;
        try {
            metamixture_core::Params p = best.p;
            // Use the largest-weight component as the reference logit.
            int ref = static_cast<int>(std::max_element(p.pi.begin(), p.pi.end()) - p.pi.begin());
            std::swap(p.mu[ref], p.mu[K-1]); std::swap(p.tau[ref], p.tau[K-1]);
            std::swap(p.pi[ref], p.pi[K-1]);
            std::vector<double> start_theta = constraints.pack(p);
            NumericVector par = Rcpp::wrap(start_theta);
            NumericVector lower(constraints.q, R_NegInf), upper(constraints.q, R_PosInf);
            const double start_floor = 1e-12;
            for (int j=0;j<constraints.nm;++j) { lower[j]=d.ymin;upper[j]=d.ymax; }
            for (int j=0;j<constraints.nt;++j) { lower[constraints.nm+j]=0;upper[constraints.nm+j]=d.upper; }
            for (double w:p.pi) if (w < start_floor) joint_start_adjusted=true;
            XPtr<metamixture_core::ConstrainedJoint> context(
                new metamixture_core::ConstrainedJoint(d, constraints), true);
            context->evaluate(start_theta);
            joint_start_ll = d.original_ll(-d.n * context->value);
            Environment stats = Environment::namespace_env("stats");
            Function optim = stats["optim"];
            List opt = optim(_["par"] = par,
                _["fn"] = InternalFunction(&metamixture_rcpp::constrained_fn),
                _["gr"] = InternalFunction(&metamixture_rcpp::constrained_gr),
                _["method"] = "L-BFGS-B", _["lower"] = lower, _["upper"] = upper,
                _["control"] = List::create(_["maxit"] = joint_maxit,
                    _["factr"] = joint_factr, _["pgtol"] = joint_pgtol),
                _["context"] = context);
            joint_code = as<int>(opt["convergence"]);
            if (opt.containsElementNamed("message") && !Rf_isNull(opt["message"]))
                joint_message = as<std::string>(opt["message"]);
            else joint_message = (joint_code == 0) ? "Converged" : "Optimizer did not report convergence";
            joint_counts = as<IntegerVector>(opt["counts"]);
            std::vector<double> theta = as<std::vector<double> >(opt["par"]), lp, rr;
            context->evaluate(theta);
            metamixture_core::Params candidate = context->unpack(theta, lp);
            // Check the exact representation that will be returned, rather
            // than relying only on optim's reported value or convergence flag.
            double candidate_ll = metamixture_core::estep(d, candidate, rr);
            joint_ll = d.original_ll(candidate_ll);
            pg_norm = 0;
            for (int j = 0; j < constraints.q; ++j) {
                double g = context->gradient[j];
                double projected = theta[j] - metamixture_core::clamp(theta[j] - g, lower[j], upper[j]);
                pg_norm = std::max(pg_norm, std::fabs(projected));
            }
            if (std::isfinite(candidate_ll) && candidate_ll >= best.ll) {
                final_p = candidate; joint_accepted = true;
            } else joint_message += "; candidate rejected because likelihood did not improve";
        } catch (const Rcpp::internal::InterruptedException&) {
            throw; // Do not suppress a user interrupt.
        } catch (const std::exception& e) {
            // A failed optional refinement must not destroy a valid ECM fit.
            joint_message = std::string("Joint refinement failed; retained ECM fit: ") + e.what();
        }
    }
    if (constraints.equal_mu) {
        std::vector<int> ix(K); std::iota(ix.begin(),ix.end(),0);
        std::stable_sort(ix.begin(),ix.end(),[&](int a,int b){return final_p.tau[a]<final_p.tau[b];});
        metamixture_core::Params tmp(K);
        for(int k=0;k<K;++k){tmp.mu[k]=final_p.mu[ix[k]];tmp.tau[k]=final_p.tau[ix[k]];tmp.pi[k]=final_p.pi[ix[k]];}
        final_p=tmp;
    } else final_p = metamixture_core::sorted(final_p);
    std::vector<double> rr;
    double final_ll_std = metamixture_core::estep(d, final_p, rr);
    double ll = d.original_ll(final_ll_std);
    NumericMatrix posterior(d.n, K);
    NumericVector effective_n(K);
    IntegerVector classification(d.n);
    for (int i = 0; i < d.n; ++i) {
        int argmax = 0;
        for (int k = 0; k < K; ++k) {
            double r = rr[static_cast<std::size_t>(i)*K+k];
            posterior(i,k) = r; effective_n[k] += r;
            if (r > rr[static_cast<std::size_t>(i)*K+argmax]) argmax = k;
        }
        classification[i] = argmax + 1;
    }
    CharacterVector notes;
    if (!best.converged) notes.push_back("Best ECM start did not meet its stopping criterion.");
    for (int k = 0; k < K; ++k) {
        if (effective_n[k] < 1.0) { notes.push_back("At least one component has effective membership < 1 study."); break; }
    }
    for (double w : final_p.pi) {
        if (w <= 1e-8) { notes.push_back("A mixing proportion is near zero; associated parameters may be unidentified."); break; }
    }
    if (d.upper == 0 && K > 1) notes.push_back("All yi are equal; the K-component representation is not identifiable.");
    if (joint_attempted && joint_code != 0) notes.push_back("Joint optimization did not report convergence; inspect $joint.");
    List estimates = metamixture_rcpp::params_list(d, final_p);
    List ecm_estimates = metamixture_rcpp::params_list(d, metamixture_core::sorted(best.p));
    NumericVector trace(best.trace.size());
    for (std::size_t j = 0; j < best.trace.size(); ++j) trace[j] = d.original_ll(best.trace[j]);
    int npar = constraints.q;
    bool converged = joint_accepted ? (joint_code == 0) : best.converged;
    Rcpp::List out = List::create(
      _["mu"] = estimates["mu"], _["tau2"] = estimates["tau2"], _["pi"] = estimates["pi"],
      _["mu_structure"] = mu, _["tau2_structure"] = tau2,
                        _["posterior"] = posterior, _["classification"] = classification,
                        _["effective_n"] = effective_n, _["logLik"] = ll,
                        _["AIC"] = -2*ll + 2*npar, _["BIC"] = -2*ll + std::log(static_cast<double>(d.n))*npar,
                        _["n"] = d.n, _["K"] = K, _["npar"] = npar, _["converged"] = converged,
                        _["ecm"] = List::create(_["mu"] = ecm_estimates["mu"], _["tau2"] = ecm_estimates["tau2"],
                                             _["pi"] = ecm_estimates["pi"], _["logLik"] = d.original_ll(best.ll),
                                               _["iterations"] = best.iterations, _["converged"] = best.converged,
                                               _["status"] = best.status, _["best_start"] = best_start, _["trace"] = trace),
                                                 _["joint"] = List::create(_["attempted"] = joint_attempted, _["accepted"] = joint_accepted,
                                                                      _["convergence"] = joint_code, _["message"] = joint_message,
                                                                      _["start_logLik"] = joint_start_ll, _["logLik"] = joint_ll,
                                                                      _["counts"] = joint_counts, _["projected_gradient_norm"] = pg_norm,
                                                                      _["starting_weights_adjusted"] = joint_start_adjusted),
                                                                      _["starts"] = DataFrame::create(_["start"] = start_number, _["logLik"] = start_ll,
                                                                                                _["iterations"] = start_it, _["converged"] = start_conv,
                                                                                                _["active_components"] = active, _["status"] = start_status),
                                                                                                _["notes"] = notes);

    if (compute_se) {
        try {
            Rcpp::List se_result = metamixture_rcpp::standard_errors(
                d, final_p, Rcpp::as<std::vector<double> >(yi),
                se_rel_step, se_boundary_tol, converged, constraints);
            out["se_mu"] = se_result["se_mu"];
            out["se_tau2"] = se_result["se_tau2"];
            out["se_pi"] = se_result["se_pi"];
            out["se"] = se_result["se"];
            out["vcov"] = se_result["vcov"];
            out["information"] = se_result["information"];
            out["se_details"] = se_result;
            out["se_status"] = converged ? "ok" : "fit_not_converged";
        } catch (const std::exception& e) {
            out["se_mu"] = Rcpp::NumericVector(K, NA_REAL);
            out["se_tau2"] = Rcpp::NumericVector(K, NA_REAL);
            out["se_pi"] = Rcpp::NumericVector(K, NA_REAL);
            out["se_status"] = std::string(e.what());
        }
    } else {
        out["se_status"] = "disabled";
    }

    out.attr("class") = Rcpp::CharacterVector::create("remix_cpp", "list");

    return out;
}
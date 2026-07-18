# Offline symbolic simplification of a tree via SymPy.
#
# `simplify(node)` walks the tree, builds an equivalent SymPy expression, runs
# `sympy.simplify`, and returns a string. It is meant for analysis of evolved
# expressions, not for anything on the hot path.
#
# `simplify(node, inputs)` additionally derives per-variable sign assumptions
# from the data (positive / nonnegative / real), unlocking stronger
# simplifications (e.g. `sqrt(x^2) -> x` when `x > 0`).

# Walk context: the SymPy handle plus the per-variable assumption lookup.
struct _Ctx
    sympy::Py
    assume::Function        # i::Int -> :positive | :nonnegative | :real
end

# --- constant / symbol construction ------------------------------------------

# Classify a variable's data into the strongest honest sign assumption.
function _assumption(vals::AbstractVector{<:Real})
    all(>(0), vals)  && return :positive
    all(>=(0), vals) && return :nonnegative
    return :real
end

function _symbol(sympy, i::Int, assume::Symbol)
    assume === :positive    && return sympy.Symbol("x$i", positive=true)
    assume === :nonnegative && return sympy.Symbol("x$i", nonnegative=true)
    return sympy.Symbol("x$i", real=true)
end

# Inject a raw numeric value as a SymPy number. Integer-valued floats become
# `Integer` (so `3.0` prints as `3`), everything else stays a `Float`.
#
# Only values up to Float64's exact-integer limit (`_NUM_CAP == 2^53`) go to
# `Integer`. Above it an integer-valued float's low digits are binary artifacts
# — e.g. `-2.9e105` expands to a 106-digit "exact" integer whose trailing digits
# are meaningless — so a compact `Float` is both honest and readable, and avoids
# the `Int64` overflow that `Int(v)` would hit. (A region-aware pass could pick
# the shortest `m*10^n` value inside the equivalence interval instead of the raw
# float; deferred — see _nice_from_region.)
function _num(sympy, v::Real)
    (isinteger(v) && abs(v) <= _NUM_CAP) && return sympy.Integer(Int(v))
    return sympy.Float(v)
end

# --- simplest constant within an equivalence region --------------------------

const _RATIONAL_DEN_CAP = 10_000    # bare rationals: denominator ceiling
const _CONST_CAP = 20               # rational-multiple forms: |num|,den ceiling
const _NUM_CAP = 2^53               # largest numerator kept exact as Float64/Int
const _CONST_PENALTY = 2            # cost of introducing a named constant

# Named irrational constants tried as `r * c`: (numeric value, SymPy expr).
# List order is the tie-break rank (earlier = preferred).
_named_constants(sympy) = [
    (Float64(π),          sympy.pi),
    (exp(1.0),            sympy.E),
    (sqrt(2.0),           sympy.sqrt(2)),
    (sqrt(3.0),           sympy.sqrt(3)),
    (sqrt(5.0),           sympy.sqrt(5)),
    ((1 + sqrt(5.0)) / 2, sympy.GoldenRatio),
]

# Smallest-denominator rational p/q with lo <= p/q <= hi, or `nothing` if none
# has denominator <= `den_cap`. Continued-fraction (Stern-Brocot) search: runs
# in O(number of CF terms), not O(value), so wide/large intervals resolve in a
# handful of steps instead of walking mediant-by-mediant. Handles sign and zero.
function _simplest_rational(lo::Real, hi::Real; den_cap::Int)
    lo <= 0 <= hi && return (0, 1)
    if hi < 0
        r = _simplest_rational(-hi, -lo; den_cap=den_cap)
        return r === nothing ? nothing : (-r[1], r[2])
    end
    # 0 < lo <= hi
    return _cf_simplest(float(lo), float(hi), den_cap, 0, 1)
end

# Simplest p/q in [lo, hi] for 0 < lo <= hi, with q <= den_cap, else `nothing`.
#
# `qb_prev`/`qb_cur` are a running Fibonacci lower bound on the denominator of
# any rational reachable at or below this level (each CF term is >= 1, so the
# convergent denominators grow at least Fibonacci-fast). Pruning the descent as
# soon as that bound passes `den_cap` bounds the depth to ~log_phi(den_cap): a
# degenerate interval (lo == hi) otherwise recurses on its own continued
# fraction forever — floating-point rounding keeps it from ever landing on an
# exact integer — and blows the stack. The bound never exceeds the true
# denominator at the same depth, so no valid result is pruned.
function _cf_simplest(lo::Float64, hi::Float64, den_cap::Int, qb_prev::Int, qb_cur::Int)
    qb_cur > den_cap && return nothing
    # The smallest-|value| integer in the interval is the simplest possible
    # rational (denominator 1). ceil(lo) is that integer when it fits under hi.
    n = ceil(lo)
    if n <= hi
        # Guard against numerators too large to hold exactly as an Int: a region
        # this wide pins no "nice" constant, so report none rather than overflow.
        (isfinite(n) && n <= _NUM_CAP) || return nothing
        return (Int(n), 1)
    end
    # No integer strictly inside, so floor(lo) == floor(hi). Peel it off and
    # recurse on the reciprocal of the fractional part (the next CF term).
    fl = floor(lo)
    sub = _cf_simplest(1 / (hi - fl), 1 / (lo - fl), den_cap, qb_cur, qb_prev + qb_cur)
    sub === nothing && return nothing
    n2, d2 = sub                       # p/q = fl + d2/n2 = (fl*n2 + d2)/n2
    n2 > den_cap && return nothing     # denominators only grow with depth
    return (Int(fl) * n2 + d2, n2)
end

# Digit-count cost of a rational's written form.
_digits(n::Integer) = ndigits(abs(n))

# Push a `r*c` (multiply=true) or `r/c` (multiply=false) candidate for the
# rational `r`, skipping ones that are absent, zero, or over the size cap.
function _push_const!(out, sympy, rank::Int, r, cexpr, multiply::Bool)
    r === nothing && return
    p, q = r
    (p == 0 || abs(p) > _CONST_CAP) && return          # 0 is just the bare case
    coef = sympy.Rational(p, q)
    expr = multiply ? coef * cexpr : coef / cexpr
    push!(out, (_digits(p) + _digits(q) + _CONST_PENALTY, rank, expr))
end

# Candidate constants for one interval [lo, hi], each as `(cost, rank, expr)`.
# Lower `(cost, rank)` wins; rank 0 (bare rational) breaks ties in its favor so
# a named constant is only chosen when it is genuinely simpler.
function _candidates(sympy, lo::Real, hi::Real, consts)
    out = Tuple{Int,Int,Py}[]

    # rank 0: bare rational p/q in the interval.
    r = _simplest_rational(lo, hi; den_cap=_RATIONAL_DEN_CAP)
    if r !== nothing
        p, q = r
        push!(out, (_digits(p) + _digits(q), 0, sympy.Rational(p, q)))
    end

    # rank k>=1: rational multiple r*c (r in [lo/cv, hi/cv]) and reciprocal
    # form s/c (s in [lo*cv, hi*cv], catching values like 1/pi).
    for (rank, (cv, cexpr)) in enumerate(consts)
        mult = _simplest_rational(lo / cv, hi / cv; den_cap=_CONST_CAP)
        _push_const!(out, sympy, rank, mult, cexpr, true)
        recip = _simplest_rational(lo * cv, hi * cv; den_cap=_CONST_CAP)
        _push_const!(out, sympy, rank, recip, cexpr, false)
    end

    return out
end

# The simplest constant lying anywhere in `regions`, as a SymPy expr, or
# `nothing` if nothing simple enough fits.
function _nice_from_region(sympy, regions::Vector{IntervalType})
    consts = _named_constants(sympy)
    best = nothing        # (cost, rank, expr)
    for iv in regions
        for cand in _candidates(sympy, iv.lo, iv.hi, consts)
            if best === nothing || (cand[1], cand[2]) < (best[1], best[2])
                best = cand
            end
        end
    end
    return best === nothing ? nothing : best[3]
end

# --- node -> SymPy -----------------------------------------------------------

# The SymPy value contributed by a node's `AddedValue`. When the value carries an
# equivalence region, prefer the simplest constant in that region over the raw
# float; otherwise fall back to the honest float.
function _offset_expr(ctx::_Ctx, add::AddedValue)
    if !isempty(add.allowed_intervals)
        nice = _nice_from_region(ctx.sympy, add.allowed_intervals)
        nice === nothing || return nice
    end
    return _num(ctx.sympy, add.value)
end

# `_core` is the node's expression WITHOUT its additive offset; `_to_sympy` then
# adds the offset uniformly. A `Constant`'s core is 0 — its whole value lives in
# the `AddedValue`, so it flows through the same offset path.
_core(ctx::_Ctx, node::Var) = _symbol(ctx.sympy, node.index, ctx.assume(node.index))
_core(ctx::_Ctx, node::Constant) = ctx.sympy.Integer(0)

function _core(ctx::_Ctx, node::BinaryNode)
    l = _to_sympy(ctx, node.left)
    r = _to_sympy(ctx, node.right)
    f = node.fun
    f === (+) && return l + r
    f === (-) && return l - r
    f === (*) && return l * r
    f === (/) && return l / r
    error("simplify: unsupported binary function $(f)")
end

function _core(ctx::_Ctx, node::UnaryNode)
    c = _to_sympy(ctx, node.child)
    sympy = ctx.sympy
    f = node.fun
    f === (-)  && return -c
    f === inv  && return sympy.Integer(1) / c
    f === sqrt && return sympy.sqrt(c)
    f === log  && return sympy.log(c)
    f === exp  && return sympy.exp(c)
    f === sin  && return sympy.sin(c)
    f === cos  && return sympy.cos(c)
    error("simplify: unsupported unary function $(f)")
end

function _to_sympy(ctx::_Ctx, node::Node)
    expr = _core(ctx, node)
    if node.addition.value != 0
        expr = expr + _offset_expr(ctx, node.addition)
    end
    return expr
end

# --- public API --------------------------------------------------------------

function _render(sympy, expr, form::Symbol)
    form === :str    && return pyconvert(String, pystr(expr))
    form === :pretty && return pyconvert(String, sympy.pretty(expr))
    form === :latex  && return pyconvert(String, sympy.latex(expr))
    error("simplify: unknown form $(form) (expected :str, :pretty, or :latex)")
end

function _simplify(node::Node, assume::Function, form::Symbol)
    sympy = pyimport("sympy")
    ctx = _Ctx(sympy, assume)
    expr = sympy.simplify(_to_sympy(ctx, node))
    return _render(sympy, expr, form)
end

simplify(node::Node; form::Symbol=:str) = _simplify(node, _ -> :real, form)

function simplify(node::Node, inputs::Vector{Vector{Float64}}; form::Symbol=:str)
    assume(i::Int) = (1 <= i <= length(inputs)) ? _assumption(inputs[i]) : :real
    return _simplify(node, assume, form)
end

# Convenience: simplify a whole Tree by delegating to its root.
simplify(tree::Tree; kwargs...) = simplify(tree.root; kwargs...)
simplify(tree::Tree, inputs::Vector{Vector{Float64}}; kwargs...) =
    simplify(tree.root, inputs; kwargs...)

export simplify

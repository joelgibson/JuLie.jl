# Coxeter Groups

A *Coxeter system* $(W, S)$ is a group $W$ together with a set of generators $S$, such that $W$ has the presentation

$W = \langle S \mid (st)^{m_{st}} = 1 \text{ for all } s, t \in S \rangle,$

where $m_{ss} = 1$ and $m_{st} \in \{0, 2, 3, 4, \ldots\}$ for $s \neq t$, with $m_{st} = 0$ meaning that there is no relation (the symbol $m_{st} = \infty$ is used in most literature on Coxeter systems). A *Coxeter group* is a group $W$ admitting a set of generators $S \subseteq W$ such that $(W, S)$ is a Coxeter system. The *rank* of a Coxeter system is $|S|$.

In JuLie, a Coxeter group will always be represented as a Coxeter system, with a finite set $S$ of generators. There is a common [abstract interface](@ref coxeter_abstract_interface) for working with Coxeter groups, and multiple different concrete implementations:

- The [Matrix representation](@ref) represents elements as pairs of matrices (a group element and its inverse) acting on the basis of simple roots. It is a baseline implementation, slow but simple, and used to double-check the other implementations.



## [Abstract interface] (@id coxeter_abstract_interface)

There are two abstract supertypes:

```@docs
CoxGrp
CoxElt
```

A minimal implementation of a Coxeter group is a group object subtyping `CoxGrp` and an element object subtyping `CoxElt`, implementing at least the following interfaces:

| Required methods                  | Brief description                                      |
|:----------------------------------|:-------------------------------------------------------|
| `rank(W::CoxGrp)`                 | The rank of the Coxeter group, a nonnegative integer   |
| `gens(W::CoxGrp)`                 | A list of length `rank(grp)`, the standard generators. |
| `one(W::CoxGrp)`                  | The identity element.                                  |
| `parent(x::CoxElt)`               | Return the group object                                |
| `left_descent_set(x::CoxElt)`     | A sub-`BitSet` of `1:rank(grp)`                        |
| `right_descent_set(x::CoxElt)`    | A sub-`BitSet` of `1:rank(grp)`                        |
| `mult_gen(s::Integer, x::CoxElt)` | Left multiplication by the generator with index `s`    |
| `mult_gen(x::CoxElt, s::Integer)` | Right multiplication by the generator with index `s`   |

In addition, the element object should provide implementations for `isequal`, `==`, and `hash`. The following methods have default implementations based on the table above:

```@docs
length(::CoxElt)
isone(::CoxElt)
lexword(::CoxElt)
*(::CoxElt, ::CoxElt)
```


## Coxeter and Cartan matrices

A Coxeter system $(W, S)$ is determined by its *Coxeter matrix* $[m_{st}]_{s, t \in S}$, a symmetric matrix with entries in $\mathbb{N}$, and $m_{ss} = 1$ for all $s \in S$. If there is no relation between the generators $s, t$, then we use the convention $m_{st} = 0$ rather than $m_{st} = \infty$.

A *Cartan matrix* is also determines a Coxeter system, however there is strictly more information in a Cartan matrix than a Coxeter matrix (many Cartan matrices determine the same Coxeter matrix). In this context, a Cartan matrix $A = [a_{st}]_{s, t \in S}$ is a matrix with real entries such that

1. Diagonal entries are $2$, and off-diagonal entries are zero or negative.
2. Off-diagonal entries satisfy $a_{st} = 0 \iff a_{ts} = 0$.
3. If $0 \leq n_{st} = a_{st} a_{ts} < 4$, then $n_{st} = 4 \cos^2(\pi / m_{st})$ for some integer $m_{st} \geq 2$.

A Cartan matrix encodes the pairings $a_{st} = \langle \alpha_s^\vee, \alpha_t \rangle$ between the simple roots and coroots in a realisation $(V^*, V)$ of the Coxeter group.


## Matrix representation

Given a Cartan matrix $A = [a_{st}]_{s, t \in S}$, let $V$ be the real vector space

```@docs
coxeter_group_mat
```

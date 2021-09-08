# Coxeter Groups

JuLie provides an [abstract interface](@ref abstract-interface) for working with Coxeter group and Coxeter group elements, and then various concrete implementations of Coxeter groups:

- The [Matrix representation](@ref) represents elements as pairs of matrices (a group element and its inverse) acting on the basis of simple roots. It is a baseline implementation, slow but simple, and used to double-check the other implementations.


## Notation

A *Coxeter system* $(W, S)$ is a group $W$ together with a set of generating involutions $S$, such that the group has presentation

$W = \langle S \mid s^2 = 1, (st)^{m_{st}} = 1 \text{ for all } s, t \in S \text{ with } m_{st} \neq 0 \rangle,$

for some integers $m_{st} = m_{ts} \in \{0, 2, 3, 4, \ldots\}$. (In the literature it is standard to use $m_{st} = \infty$ to mean no relation between $s$ and $t$, wheras in JuLie we use $m_{st} = 0$). In JuLie a Coxeter group will always be represented by a Coxeter system, with the generators numbered $\{1, \ldots, r\}$ where $r$ is the *rank* of the Coxeter system.

A *Coxeter matrix* $[m_{st}]_{s, t \in S}$ is a square symmetric matrix such that $m_{ss} = 1$ for all $s \in S$, and $m_{st} \in \{0, 2, 3, 4, \ldots\}$ for all $s \neq t$. Each Coxeter matrix determines a unique Coxeter system. A Coxeter matrix or Coxeter group is called *crystallographic* if $m_{st} \in \{0, 2, 3, 4, 6\}$ for all $s \neq t$.

A *Cartan matrix* is a square matrix $A = [a_{st}]_{s, t \in S}$ with real entries such that $a_{ss} = 2$ for all $s \in S$, $a_{st} \leq 0$ for all $s \neq t$, $a_{st} = 0 \iff a_{ts} = 0$, and if $a_{st}a_{ts} < 4$ for $s \neq t$, then $a_{st}a_{ts} = 4 \cos^2(\pi/m_{st})$ for some integers $m_{st} \geq 2$. A Cartan matrix determines a Coxeter matrix via the rules $m_{st} = 0$ if $a_{st}a_{ts} \geq 4$, and otherwise $4 \cos^2(\pi/m_{st}) = a_{st} a_{ts}$. Note that $a_{st}a_{ts} = 0, 1, 2, 3$ implies $m_{st} = 2, 3, 4, 6$ respectively. A Cartan matrix is called *crystallographic* if it has integer entries, and in this case the resulting Coxeter system is also crystallographic, but a crystallographic Coxeter group may admit non-crystallographic Cartan matrices.

A Cartan matrix records some essential information about a faithful reflection representation of a Coxeter group. Let $V$ be a finite-dimensional real vector space and suppose there exist sets $\{\alpha_s\}_{s \in S} \subseteq V$ and $\{\alpha_s^\vee\}_{s \in S} \subseteq V^*$ such that $\langle \alpha_s^\vee, \alpha_t \rangle = a_{st}$. The reflections $r_s \colon V \to V$ defined by $r_s(v) = v - \langle \alpha_s^\vee, v \rangle \alpha_s$ then form a representation of $W$ inside $\operatorname{GL}(V)$, and provided that either the *simple roots* $\{\alpha_s\}$ or the *simple coroots* $\{\alpha_t\}$ are linearly independent, this representation is faithful.

Taking $V$ to be the vector space with basis $\{\alpha_s\}_{s \in S}$ always yields such a faithful representation, for any Cartan matrix. Finally, given any Coxeter matrix $[m_{st}]_{s, t \in S}$, a compatible Cartan matrix is made by setting $a_{st} = -2 \cos(\pi / m_{st})$ for all $m_{st} \geq 1$, and $a_{st} = 2$ if $m_{st} = 0$: this resulting representation $V$ is called the *geometric representation* of the Coxeter system.


## [Abstract interface](@id abstract-interface)

There are two abstract supertypes, one to represent Coxeter groups (i.e. Coxeter systems), and another to represent group elements:

```@docs
CoxGrp
CoxElt
```

A minimal implementation of a Coxeter group is a group object subtyping `CoxGrp` and an element object subtyping `CoxElt`, implementing at least the following interfaces:

| Required methods                  | Brief description                                      |
|:----------------------------------|:-------------------------------------------------------|
| `rank(W::CoxGrp)`                 | The rank of the Coxeter group, a nonnegative integer.  |
| `gens(W::CoxGrp)`                 | A list of length `rank(grp)`, the standard generators. |
| `one(W::CoxGrp)`                  | The identity element.                                  |
| `parent(x::CoxElt)`               | Return the group object.                               |
| `left_descent_set(x::CoxElt)`     | A sub-`BitSet` of `1:rank(grp)`                        |
| `right_descent_set(x::CoxElt)`    | A sub-`BitSet` of `1:rank(grp)`                        |
| `mult_gen(s::Integer, x::CoxElt)` | Left multiplication by the generator with index `s`    |
| `mult_gen(x::CoxElt, s::Integer)` | Right multiplication by the generator with index `s`   |

In addition, the element object should provide implementations for `isequal`, `==`, and `hash`. The following methods have default implementations based on the table above:

```@docs
length(::CoxElt)
isone(::CoxElt)
iterate(::CoxElt)
lexword(::CoxElt)
*(::CoxElt, ::CoxElt)
^(::CoxElt, ::Integer)
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

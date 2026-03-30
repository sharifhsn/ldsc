// sn-theorems.typ — Theorem environments matching sn-jnl.cls styles
//
// Three theorem styles from sn-jnl.cls:
//   thmstyleone:   bold head, italic body  (theorem, proposition)
//   thmstyletwo:   italic head, normal body (example, remark)
//   thmstylethree: bold head, normal body  (definition)
// Plus proof environment with QED square.
//
// Usage:
//   #import "sn-theorems.typ": theorem, proposition, example, remark, definition, proof

// ---------------------------------------------------------------------------
// Shared counter per top-level heading
// ---------------------------------------------------------------------------
#let _thm-counter = counter("sn-theorem")

#let _make-theorem(supplement, head-style, body-style) = {
  (title: none, body) => {
    _thm-counter.step()
    v(18pt)
    block(width: 100%, {
      set text(size: 9pt)
      // Head
      {
        set text(..head-style)
        [#supplement ]
        context _thm-counter.display("1")
        if title != none [ (#title)]
        [. ]
      }
      // Body
      {
        set text(..body-style)
        body
      }
    })
    v(18pt)
  }
}

// thmstyleone: bold head, italic body
#let theorem = _make-theorem(
  "Theorem",
  (weight: "bold"),
  (style: "italic"),
)

#let proposition = _make-theorem(
  "Proposition",
  (weight: "bold"),
  (style: "italic"),
)

// thmstyletwo: italic head, normal body
#let example = _make-theorem(
  "Example",
  (style: "italic"),
  (:),
)

#let remark = _make-theorem(
  "Remark",
  (style: "italic"),
  (:),
)

// thmstylethree: bold head, normal body
#let definition = _make-theorem(
  "Definition",
  (weight: "bold"),
  (:),
)

// ---------------------------------------------------------------------------
// Proof environment with QED square
// ---------------------------------------------------------------------------
#let proof(title: none, body) = {
  v(18pt)
  block(width: 100%, {
    set text(size: 9pt)
    {
      set text(style: "italic")
      [Proof]
      if title != none [ (#title)]
      [. ]
    }
    body
    h(1fr)
    sym.square.stroked
  })
  v(18pt)
}

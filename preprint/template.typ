// bioRxiv-style preprint template
// Based on preprintx layout with bug fixes and font fallbacks

#let preprint(
  title: none,
  authors: (),  // list of (name: str, affiliations: str, corresponding: bool)
  affils: (:),  // "1": "Institution Name"
  abstract: none,
  keywords: (),
  correspondence: none,
  date: none,
  doc,
) = {
  // Font stack with fallbacks
  let sans = ("Liberation Sans",)
  let serif = ("STIX Two Text", "Linux Libertine", "Libertinus Serif", "New Computer Modern")

  set text(font: serif, size: 10pt)
  show heading: set text(font: sans)
  show heading.where(level: 1): it => {
    set text(11pt, weight: "bold")
    pad(bottom: 8pt, top: 4pt, it)
  }
  show heading.where(level: 2): it => {
    set text(9pt, weight: "bold")
    box(it) + "."
  }
  show heading.where(level: 3): it => {
    set text(9pt)
    emph(box(it) + ".")
  }

  show bibliography: set text(font: sans, size: 8pt)
  show figure.caption: set text(font: sans, size: 8pt)
  show figure: it => {
    place(auto, float: true, scope: "parent")[
      #it.body
      #it.caption
    ]
  }

  // Footer with author names
  let footer-name = if authors.len() == 1 {
    authors.at(0).name
  } else if authors.len() == 2 {
    authors.at(0).name.split(",").at(0) + " & " + authors.at(1).name.split(",").at(0)
  } else {
    authors.at(0).name.split(",").at(0) + " et al."
  }

  set page(
    paper: "us-letter",
    margin: 0.75in,
    footer: context align(right + horizon)[
      #text(font: sans, size: 7pt)[
        #footer-name #sym.bar.v _Preprint_ #sym.bar.v #counter(page).display("1")
      ]
    ],
    columns: 2,
  )

  // Title block (spans both columns)
  place(top + center, float: true, scope: "parent")[
    #set align(center)
    #text(24pt, weight: "bold", title)

    #v(8pt)
    #text(font: sans, size: 8pt, weight: "bold")[
      #authors.map(a => {
        let parts = a.name.split(",").map(s => s.trim())
        let display = if parts.len() >= 2 { parts.at(1) + " " + parts.at(0) } else { a.name }
        [#display#super(a.affiliations)]
      }).join(", ", last: " and ")
    ]

    #v(4pt)
    #text(font: sans, size: 7pt)[
      #affils.pairs().map(p => [#super(p.at(0))#p.at(1)]).join([\ ])
    ]

    #if date != none {
      v(4pt)
      text(font: sans, size: 7pt, style: "italic", date)
    }
  ]

  // Abstract
  set align(left)
  set par(justify: true)

  text(font: serif, size: 9pt, weight: "bold", abstract)

  v(4pt)
  text(font: sans, size: 7pt, weight: "bold")[
    #keywords.join([ #sym.bar.v ])
  ]

  if correspondence != none {
    linebreak()
    text(font: sans, size: 7pt, style: "italic", weight: "bold")[\u{2709} #correspondence]
  }

  v(8pt)

  doc
}

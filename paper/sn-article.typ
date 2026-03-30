// sn-article.typ — Springer Nature journal article template for Typst
// Replicates the sn-jnl LaTeX class (v2.1, April 2023)
//
// Usage:
//   #import "sn-article.typ": sn-article, wide-figure, wide-table, appendix, backmatter-heading
//   #show: sn-article.with(title: [...], authors: (...), ...)

// ---------------------------------------------------------------------------
// Bibliography style mapping: sn-jnl option → CSL style name
//
// These map to Typst's built-in CSL styles, which come from the official CSL
// repository. They are close approximations of the corresponding .bst files
// but may differ in minor formatting details. Springer reformats accepted
// papers during production, so these are sufficient for submission.
// For exact matching, replace the string value with a path to a custom .csl file.
// ---------------------------------------------------------------------------
#let _bib-style-map = (
  "sn-basic": "springer-basic-author-date",
  "sn-nature": "nature",
  "sn-mathphys": "springer-mathphys",
  "sn-aps": "american-physics-society",
  "sn-vancouver": "vancouver",
  "sn-apa": "apa",
  "sn-chicago": "chicago-author-date",
)

// ---------------------------------------------------------------------------
// Helper: format author running head (abbreviated for > 2 authors)
// ---------------------------------------------------------------------------
#let _author-running-head(authors) = {
  let names = authors.map(a => {
    if a.at("surname", default: none) != none { a.surname }
    else {
      let parts = a.name.split(" ")
      parts.last()
    }
  })
  if names.len() == 0 { return [] }
  if names.len() == 1 { return names.first() }
  if names.len() == 2 { return names.at(0) + " and " + names.at(1) }
  return names.first() + " et al."
}

// ---------------------------------------------------------------------------
// Helper: render ORCID icon (inline green icon linking to profile)
// ---------------------------------------------------------------------------
#let _orcid-icon(id) = {
  link(
    "https://orcid.org/" + id,
    box(
      height: 0.9em,
      baseline: 15%,
      image("fig/orcid.svg"),
    ),
  )
}

// ---------------------------------------------------------------------------
// Helper: render a single author with affiliation superscripts
// ---------------------------------------------------------------------------
#let _render-author(author, show-affil-nums: true, total-authors: 1) = {
  let name = if author.at("given-name", default: none) != none and author.at("surname", default: none) != none {
    author.given-name + " " + author.surname
  } else {
    author.name
  }

  name

  // ORCID icon
  if author.at("orcid", default: none) != none {
    h(2pt)
    _orcid-icon(author.orcid)
  }

  // Affiliation superscripts
  if show-affil-nums and author.at("affiliations", default: ()).len() > 0 {
    super(author.affiliations.join(","))
  }

  // Equal contribution dagger
  if author.at("equal-contribution", default: false) {
    super(sym.dagger)
  }

  // Corresponding author asterisk
  if author.at("corresponding", default: false) and total-authors > 1 {
    super("*")
  }
}

// ---------------------------------------------------------------------------
// Wide figure / wide table — span both columns in 2-column mode
// ---------------------------------------------------------------------------
#let wide-figure(caption: none, placement: auto, label: none, ..args) = {
  figure(
    ..args,
    caption: caption,
    placement: placement,
    scope: "parent",
  )
}

#let wide-table(caption: none, placement: auto, label: none, ..args) = {
  figure(
    ..args,
    caption: caption,
    placement: placement,
    scope: "parent",
    kind: table,
  )
}

// ---------------------------------------------------------------------------
// Appendix: resets counters, switches numbering to A.1
// ---------------------------------------------------------------------------
#let appendix(body) = {
  counter(heading).update(0)
  set heading(numbering: "A.1")
  // Reset figure and table counters
  counter(figure.where(kind: image)).update(0)
  counter(figure.where(kind: table)).update(0)
  body
}

// ---------------------------------------------------------------------------
// Back-matter heading: unnumbered bold heading (e.g. Acknowledgments)
// ---------------------------------------------------------------------------
#let backmatter-heading(title) = {
  heading(level: 1, numbering: none, title)
}

// ---------------------------------------------------------------------------
// Main template function
// ---------------------------------------------------------------------------
#let sn-article(
  // Metadata
  title: none,
  short-title: none,
  subtitle: none,
  article-type: none,

  // Authors & affiliations
  authors: (),
  affiliations: (),

  // Content metadata
  abstract: none,
  keywords: (),
  pacs: (),

  // Layout
  columns: 1,
  referee: false,
  line-numbers: false,

  // Bibliography
  bib-style: "sn-basic",
  bib-file: none,

  // Fonts
  font-size: 10pt,
  font-family: ("Libertinus Serif", "New Computer Modern", "PT Serif"),

  // Document body
  doc,
) = {
  // -------------------------------------------------------------------------
  // Resolve derived values
  // -------------------------------------------------------------------------
  let running-title = if short-title != none { short-title } else if title != none { title } else { [] }
  let author-head = _author-running-head(authors)
  let csl-style = _bib-style-map.at(bib-style, default: bib-style)
  let is-twocol = columns == 2

  // -------------------------------------------------------------------------
  // Page geometry — from sn-jnl.cls geometry package calls
  // Single column: text={31pc,194.25mm}, top=26mm, headheight=5.5pt, headsep=5.6mm
  // Double column: text={160mm,216mm}, top=26mm, headheight=12pt, headsep=5.15mm
  // -------------------------------------------------------------------------
  let text-width = if is-twocol { 160mm } else { 31 * 12pt }
  let top-margin = if is-twocol { 26mm + 12pt + 5.15mm } else { 26mm + 5.5pt + 5.6mm }
  let text-height = if is-twocol { 216mm } else { 194.25mm }
  let bottom-margin = 297mm - top-margin - text-height

  let _header = context {
    let pn = counter(page).get().first()
    if pn > 1 {
      set text(size: 8pt, style: "italic")
      if calc.odd(pn) {
        h(1fr)
        running-title
      } else {
        author-head
        h(1fr)
      }
      v(-3pt)
      line(length: 100%, stroke: 0.2mm)
    }
  }

  let _footer = context {
    align(center, text(size: 8pt, counter(page).display()))
  }

  set page(
    paper: "a4",
    columns: columns,
    margin: (
      top: top-margin,
      bottom: bottom-margin,
      left: (210mm - text-width) / 2 + 3mm,
      right: (210mm - text-width) / 2 - 3mm,
    ),
    header: _header,
    footer: _footer,
  )

  // -------------------------------------------------------------------------
  // Global typography
  // -------------------------------------------------------------------------
  set text(font: font-family, size: font-size, lang: "en")
  set par(
    justify: true,
    leading: if referee { 1.2em } else { 0.65em },
    first-line-indent: (amount: 1.5em, all: true),
  )

  // Line numbering (referee/review mode)
  if line-numbers {
    set par.line(numbering: "1")
  }

  // -------------------------------------------------------------------------
  // Heading styles — from sn-jnl.cls \sectionfont etc.
  // Section: 14bp bold, 12pt above, 9pt below
  // Subsection: 12bp bold, 12pt above, 6pt below
  // Subsubsection: 11bp bold, 12pt above, 6pt below
  // -------------------------------------------------------------------------
  set heading(numbering: "1.1.1")

  show heading.where(level: 1): it => {
    set text(size: 14pt, weight: "bold")
    set par(first-line-indent: 0pt)
    v(12pt)
    block(it)
    v(9pt)
  }

  show heading.where(level: 2): it => {
    set text(size: 12pt, weight: "bold")
    set par(first-line-indent: 0pt)
    v(12pt)
    block(it)
    v(6pt)
  }

  show heading.where(level: 3): it => {
    set text(size: 11pt, weight: "bold")
    set par(first-line-indent: 0pt)
    v(12pt)
    block(it)
    v(6pt)
  }

  // -------------------------------------------------------------------------
  // Caption formatting — from sn-jnl.cls: 8bp/9.5bp, bold "Fig." / "Table"
  // -------------------------------------------------------------------------
  set figure(placement: auto)
  set figure.caption(separator: [ ])
  show figure.caption: it => {
    set text(size: 8pt)
    set par(leading: 0.45em, first-line-indent: 0pt)
    [*#it.supplement #context it.counter.display(it.numbering)*#it.separator#it.body]
  }

  // Figure supplement: "Fig." per Springer style
  set figure(supplement: [Fig.])
  show figure.where(kind: table): set figure(supplement: [Table])

  // -------------------------------------------------------------------------
  // Equations: numbered
  // -------------------------------------------------------------------------
  set math.equation(numbering: "(1)")

  // -------------------------------------------------------------------------
  // Lists: compact spacing like sn-jnl.cls
  // -------------------------------------------------------------------------
  set list(indent: 1.5em, spacing: 0pt)
  set enum(indent: 1.5em, spacing: 0pt)

  // =========================================================================
  // FRONT MATTER
  // =========================================================================
  {
    set par(first-line-indent: 0pt)

    // Article type badge — 8bp on gray background
    if article-type != none {
      block(
        fill: luma(200),
        inset: (x: 4pt, y: 2pt),
        text(size: 8pt, upper(article-type)),
      )
      v(20pt)
    }

    // Title — 17bp/22.5bp per Titlefont
    if title != none {
      block(
        width: 100%,
        {
          set text(size: 17pt)
          set par(leading: 0.7em)
          title
        },
      )
    }

    // Subtitle — 14bp/16.5bp
    if subtitle != none {
      v(9pt)
      block(width: 100%, text(size: 14pt, subtitle))
    }

    // Authors — 12bp/14.5bp bold
    if authors.len() > 0 {
      v(20pt)
      let show-nums = affiliations.len() > 1
      block(width: 100%, {
        set text(size: 12pt, weight: "bold")
        for (i, author) in authors.enumerate() {
          if i > 0 {
            if i == authors.len() - 1 [~and ] else [, ]
          }
          _render-author(author, show-affil-nums: show-nums, total-authors: authors.len())
        }
      })
    }

    // Affiliations — 11bp/13.5bp
    if affiliations.len() > 0 {
      v(7pt)
      for affil in affiliations {
        block(width: 100%, {
          set text(size: 11pt)
          if affiliations.len() > 1 {
            super(affil.id)
          }
          let parts = ()
          if affil.at("department", default: none) != none {
            parts.push(affil.department)
          }
          if affil.at("institution", default: none) != none {
            parts.push(affil.institution)
          }
          if affil.at("address", default: none) != none {
            parts.push(affil.address)
          }
          parts.join(", ")
        })
      }
    }

    // Corresponding author email + equal contribution notes
    {
      let corr-authors = authors.filter(a => a.at("corresponding", default: false))
      if corr-authors.len() > 0 {
        v(24pt)
        block(width: 100%, {
          set text(size: 9pt)
          if authors.len() > 1 [\*]
          [Corresponding author(s). E-mail(s): ]
          corr-authors.map(a => a.at("email", default: "")).filter(e => e != "").join(", ")
        })
      }

      let equal-authors = authors.filter(a => a.at("equal-contribution", default: false))
      if equal-authors.len() > 0 {
        v(4pt)
        block(width: 100%, {
          set text(size: 9pt)
          [#sym.dagger These authors contributed equally to this work.]
        })
      }
    }

    // Abstract — 9bp/11bp, 24pt left/right indent
    if abstract != none {
      v(24pt)
      block(
        inset: (left: 24pt, right: 24pt),
        {
          set text(size: 9pt)
          set par(leading: 0.55em)
          align(center, text(weight: "bold", [Abstract]))
          v(3pt)
          abstract
        },
      )
    }

    // Keywords — 8bp/9.5bp
    if keywords.len() > 0 {
      v(10pt)
      block(
        inset: (left: 24pt, right: 24pt),
        {
          set text(size: 8pt)
          [*Keywords:* #keywords.join(" · ")]
        },
      )
    }

    // PACS codes
    if pacs.len() > 0 {
      v(6pt)
      block(
        inset: (left: 24pt, right: 24pt),
        {
          set text(size: 8pt)
          [*PACS:* #pacs.join(" · ")]
        },
      )
    }

    v(36pt)
  }

  // =========================================================================
  // DOCUMENT BODY
  // =========================================================================
  doc

  // =========================================================================
  // BIBLIOGRAPHY
  // =========================================================================
  if bib-file != none {
    bibliography(bib-file, style: csl-style)
  }
}

# Typst Documentation Reference for Template Building

**Last Updated:** 2026-03-16
**Crawled Pages:** 13 core documentation pages + 2 Universe packages
**Purpose:** Comprehensive reference for Typst template development — function signatures, parameters, code examples, and migration guidance.

---

## Table of Contents

1. [LaTeX to Typst Migration](#latex-to-typst-migration)
2. [Template Authoring Architecture](#template-authoring-architecture)
3. [Page Layout & Setup](#page-layout--setup)
4. [Core Functions Reference](#core-functions-reference)
5. [Styling with Set & Show Rules](#styling-with-set--show-rules)
6. [Universe Packages](#universe-packages)

---

## LaTeX to Typst Migration

### Syntax Equivalents

| Element | LaTeX | Typst |
|---------|-------|-------|
| Strong/bold | `\textbf{strong}` | `*strong*` |
| Emphasis/italic | `\emph{emphasis}` | `_emphasis_` |
| Hyperlink | `\url{url}` | `url/` (auto-detected) |
| Label | `\label{intro}` | `<intro>` |
| Reference | `\ref{intro}` | `@intro` |
| Citation | `\cite{key}` | `@key` |
| Inline code | `lstlisting` | `` `code` `` |
| Bullet list | `itemize` environment | `- Item` |
| Numbered list | `enumerate` environment | `+ Item` |
| Terminology list | `description` | `/ Term: Description` |

### Heading Syntax

Replace LaTeX `\section`, `\subsection` with equals-based hierarchy:
```typst
= Level 1 Heading
== Level 2 Heading
=== Level 3 Heading
```

### Code vs. Markup Mode

- **Markup mode** (default): Text composition, uses markdown-like syntax
- **Code mode**: Programming logic, entered with `#` prefix
- **Content arguments**: Use square brackets `[...]` not curly braces

```typst
#text(weight: "bold")[This is bold content]
#underline([emphasized text])
```

### Common Function Equivalents

| LaTeX Package | Typst Function | Notes |
|---------------|---|---|
| `\textbf`, `\textit` | `#text(weight:)`, `#text(style:)` | Replaces `\bfseries`, `\itshape` declarations |
| `graphicx` | `#image()` | For inserting graphics |
| `tabularx` | `#table()` | For structured tables |
| `geometry`, `fancyhdr` | `#page()` | Single function for all page setup |
| `xcolor` | `#set text(fill: rgb(...))` | Color via fill parameter |
| `hyperref` | `#link()` | Inline hyperlinks |

### Set Rules vs. LaTeX Declarations

LaTeX:
```latex
{\bfseries This text is bold}
```

Typst:
```typst
#set text(weight: "bold")
This text is bold until scope ends.
```

Set rules modify all following content within their scope, replacing LaTeX's declarative commands like `\bfseries`.

### Template Structure

Replace `\documentclass{article}` with template functions:

```typst
#import "conf.typ": conf
#show: conf.with(title: [...], authors: [...])
```

### Bibliography

```typst
This was noted by pirates. @arrgh

Multiple sources: @arrgh @netwok

#bibliography("works.bib")
```

Supports BibLaTeX `.bib` files and native Hayagriva `.yaml`/`.yml` format.

### Achieving LaTeX Aesthetic

```typst
#set page(margin: 1.75in)
#set par(
  leading: 0.55em,
  spacing: 0.55em,
  first-line-indent: 1.8em,
  justify: true,
)
#set text(font: "New Computer Modern")
#show raw: set text(font: "New Computer Modern Mono")
#show heading: set block(above: 1.4em, below: 1em)
```

---

## Template Authoring Architecture

### Basic Template with show/set Rules

```typst
#let template(doc) = [
  #set text(font: "Inria Serif")
  #show "something cool": [Typst]
  #doc
]

#show: template
```

The document content is passed as parameter `doc` to the template function.

### Parameterized Templates

```typst
#let conf(
  authors: (),
  abstract: [],
  doc,
) = {
  set page(paper: "us-letter", columns: 2)
  set text(font: "Libertinus Serif", size: 11pt)
  doc
}

#show: conf.with(
  authors: (/* author list */),
  abstract: lorem(80),
)
```

Use `.with()` method to pre-populate named arguments without verbose closures.

### Variables & Reusability

```typst
#let ipa = [taɪpst]
#let amazed(term, color: blue) = {
  text(color, box[✨ #term ✨])
}
```

Store repeated content in variables for consistency. Functions accept positional and named parameters with defaults.

### Dynamic Layouts

```typst
#let count = authors.len()
#let ncols = calc.min(count, 3)
#grid(
  columns: (1fr,) * ncols,
  row-gutter: 24pt,
  ..authors.map(author => [
    #author.name \
    #author.affiliation \
    #link("mailto:" + author.email)
  ]),
)
```

- Use array methods (`map()`, `len()`) to process collections
- Spread operator (`..`) unpacks arrays as separate arguments
- Calculate dynamic column counts based on content

### Modular Organization

```typst
#import "conf.typ": conf
#show: conf.with(authors: [...])
```

Separate templates into dedicated files. Import and apply via `.with()` for pre-populated configuration.

---

## Page Layout & Setup

### page() Function

**Default paper size:** `"a4"`

```typst
#set page("us-letter")
#set page(width: 3cm, height: 4cm)
```

### Margin Configuration

**Accepts scalar or dictionary:**
```typst
#set page(margin: 2cm)  // All sides

#set page(margin: (
  top: 3cm,
  bottom: 2cm,
  x: 1.5cm,  // left and right
  y: 1.75cm, // top and bottom
))

#set page(margin: (
  inside: 2.5cm,  // Binding edge
  outside: 2cm,
  y: 1.75cm,
))
```

**Margin keys:**
- `top`, `bottom`, `left`, `right` — individual sides
- `x` — left and right together
- `y` — top and bottom together
- `inside`/`outside` — for binding
- `rest` — fallback for unspecified sides

### Page Numbering

**Pattern-based:**
```typst
#set page(numbering: "1")           // Arabic: 1, 2, 3
#set page(numbering: "i")           // Roman: i, ii, iii
#set page(numbering: "I")           // Roman: I, II, III
#set page(numbering: "1 of 1")      // Current of total
#set page(numbering: "— 1 —")       // Custom format
#set page(number-align: center + bottom)
```

**Counter manipulation:**
```typst
#counter(page).update(1)            // Reset to 1
#counter(page).update(n => n + 5)   // Skip 5 pages
```

### Headers & Footers

**Static content:**
```typst
#set page(header: [_Title_ #h(1fr) Organization])
```

**Context-dependent (page-aware):**
```typst
#set page(header: context {
  if counter(page).get().first() > 1 {
    [Header text]
  }
})
```

**Positioning:**
- `header-ascent`: Control top margin of header
- `footer-descent`: Control bottom margin of footer
- Headers are bottom-aligned by default
- Use `context` keyword to access page-dependent counters

### Columns

```typst
#set page(columns: 2)
#set columns(gutter: 12pt)  // Spacing between columns
#colbreak()                  // Manual column break
```

### Special Attributes

```typst
#set page(fill: rgb("444352"))        // Background color
#set page(background: [...])          // Repeat behind content
#set page(foreground: [...])          // Repeat in front
#set page(flipped: true)              // Single landscape page
#set page(binding: left)              // Binding edge
```

Typst reverts to set rule settings after single-page modifications (e.g., `flipped: true`).

---

## Core Functions Reference

### bibliography()

**Signature:**
```typst
bibliography(
  sources: str | bytes | array,
  title: none | auto | content = auto,
  full: bool = false,
  style: str | bytes = "ieee",
) -> content
```

**Parameters:**
- `sources` (required): Path(s) to `.bib` (BibLaTeX) or `.yaml`/`.yml` (Hayagriva) files
- `title`: Bibliography heading (default: language-appropriate auto)
  - `none` — no heading
  - `auto` — default language heading
  - `content` — custom heading text
- `full`: Include all sources or only cited ones (default: `false`)
- `style`: Citation style name or path (default: `"ieee"`)
  - Built-in: `"ieee"`, `"apa"`, `"chicago"`, etc.
  - Custom: path to CSL file

**Example:**
```typst
@piratesreference

#bibliography("works.bib", full: false, style: "apa")
```

---

### figure()

**Signature:**
```typst
figure(
  content,
  alt: none | str,
  placement: none | auto | alignment,
  scope: str,
  caption: none | content,
  kind: auto | str | function,
  supplement: none | auto | content | function,
  numbering: none | str | function,
  gap: length,
  outlined: bool,
) -> content
```

**Parameters:**

| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `body` | content | Required | Main content (often an image) |
| `alt` | none\|str | none | Accessibility description |
| `placement` | none\|auto\|alignment | none | Float behavior: `top`, `bottom`, `auto` |
| `scope` | str | `"column"` | Containment: `"column"` or `"parent"` |
| `caption` | none\|content | none | Figure caption text |
| `kind` | auto\|str\|function | auto | Figure type; auto-detects images, tables, code |
| `supplement` | none\|auto\|content\|function | auto | Label prefix (e.g., "Figure", "Table") |
| `numbering` | none\|str\|function | `"1"` | Numbering pattern or function |
| `gap` | length | `0.65em` | Vertical space between body and caption |
| `outlined` | bool | true | Include in table of figures |

**Key Features:**
- Auto-detects content type (image → Figure counter, table → Table counter)
- Floating figures use `placement` parameter
- `kind` parameter customizes counter and supplement

**Caption Sub-element:**
```typst
figure.caption(
  position: alignment,           // top or bottom
  separator: auto | content,
  body: content,
) -> content
```

**Example:**
```typst
#figure(
  image("diagram.png", width: 80%),
  caption: [Diagram with caption],
  placement: bottom,
  kind: "image",
)
```

---

### heading()

**Signature:**
```typst
heading(
  level: auto | int,
  depth: int,
  offset: int,
  numbering: none | str | function,
  supplement: none | auto | content | function,
  outlined: bool,
  bookmarked: auto | bool,
  hanging-indent: auto | length,
  content,
) -> content
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `level` | auto\|int | auto | Heading depth (starts at 1); auto = `offset + depth` |
| `depth` | int | 1 | Relative depth within section |
| `offset` | int | 0 | Base level offset |
| `numbering` | none\|str\|function | none | Pattern like `"I."` for Roman, `"1."` for Arabic |
| `supplement` | none\|auto\|content\|function | auto | Label prefix (e.g., "Chapter") |
| `outlined` | bool | true | Include in outline/TOC |
| `bookmarked` | auto\|bool | auto | Create PDF bookmark |
| `hanging-indent` | auto\|length | auto | Indent for wrapped lines |
| `body` | content | Required | Heading text |

**Markdown Syntax:**
```typst
= Level 1
== Level 2
=== Level 3
```

**Programmatic Creation:**
```typst
#heading(level: 1, numbering: "I.")[Heading Text]
```

**Styling:**
```typst
#show heading: set block(above: 1.4em, below: 1em)
#show heading.where(level: 1): set text(size: 20pt)
```

---

### page()

**Signature:**
```typst
page(
  paper: str,
  width: auto | length,
  height: auto | length,
  flipped: bool,
  margin: auto | relative | dictionary,
  binding: auto | alignment,
  columns: int,
  fill: none | auto | color | gradient | tiling,
  numbering: none | str | function,
  supplement: none | auto | content,
  number-align: alignment,
  header: none | auto | content,
  header-ascent: relative,
  footer: none | auto | content,
  footer-descent: relative,
  background: none | content,
  foreground: none | content,
  body: content,
) -> content
```

**Common Parameters:**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `paper` | `"a4"` | A0–A11, `"us-letter"`, `"us-legal"`, business cards, etc. |
| `width` | `595.28pt` | Settable; `auto` to fit content |
| `height` | `841.89pt` | Settable; `auto` for dynamic breaks |
| `margin` | `auto` | Dict: `{top, right, bottom, left, inside, outside, x, y, rest}` |
| `numbering` | `none` | Pattern: `"1"`, `"i"`, `"1 of 1"`, or custom function |
| `number-align` | `center` | Position of page number |
| `header` | `auto` | Fills top margin; shows number if numbering set |
| `footer` | `auto` | Fills bottom margin |
| `fill` | `auto` | Background color; `none` = transparent |

**Examples:**
```typst
#set page(width: 3cm, height: 4cm, margin: (x: 8pt, y: 4pt))
#set page(numbering: "1 of 1", number-align: center + bottom)
#set page(fill: rgb("444352"))
#set text(fill: rgb("fdfdfd"))
```

---

### text()

**Signature (abridged):**
```typst
text(
  font, fallback, style, weight, stretch, size,
  fill, stroke, tracking, spacing,
  ...other properties...
  body,
) -> content
```

**Key Parameters:**

| Parameter | Type | Default | Options |
|-----------|------|---------|---------|
| `font` | str\|array\|dict | `"libertinus serif"` | Font name, multiple fonts (fallback), dict config |
| `size` | length | `11pt` | Any length (pt, em, %, etc.) |
| `weight` | int\|str | `"regular"` | "thin", "extralight", "light", "regular", "medium", "semibold", "bold", "extrabold", "black" (or 100–900) |
| `style` | str | `"normal"` | "normal", "italic", "oblique" |
| `fill` | color\|gradient\|tiling | black | Text color |
| `lang` | str | `"en"` | Language code (affects hyphenation, etc.) |

**Examples:**
```typst
#set text(18pt)
With a set rule, all following text is 18pt.

#text(blue)[Blue text]
#text(weight: "bold", size: 14pt)[Bold, 14pt text]
#text(font: "DejaVu Sans Mono")[Monospace]
```

---

### par()

**Signature:**
```typst
par(
  leading: length,
  spacing: length,
  justify: bool,
  justification-limits: dictionary,
  linebreaks: auto | str,
  first-line-indent: length | dictionary,
  hanging-indent: length,
  body: content,
) -> content
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `leading` | length | `0.65em` | Line height (bottom of one line to top of next) |
| `spacing` | length | `1.2em` | Space between consecutive paragraphs |
| `justify` | bool | `false` | Full justification across width |
| `linebreaks` | auto\|str | `auto` | Algorithm: `"simple"` or `"optimized"` |
| `first-line-indent` | length\|dict | `(amount: 0pt, all: false)` | Initial line indent; dict: `{amount, all}` |
| `hanging-indent` | length | `0pt` | Indent all lines except the first |
| `justification-limits` | dict | `(spacing: (min: 66.67%, max: 150%), ...)` | Word/character spacing limits |
| `body` | content | Required | Paragraph text |

**Dictionary Form:**
```typst
#set par(first-line-indent: (amount: 1em, all: false))
```

**Examples:**
```typst
#set par(leading: 1.2em, spacing: 1em, justify: true)
#set par(first-line-indent: 1em)
#set par(hanging-indent: 1.5em)
```

---

### footnote()

**Signature:**
```typst
footnote(numbering, body) -> content
```

**Parameters:**
- `numbering`: Pattern (str) or function; default `"1"`
  - Settable via `#set footnote(numbering: "a")`
- `body`: Content of the footnote (required)

**Key Features:**
- Automatically creates superscript links
- Displays at page bottom
- Sequential numbering throughout doc (resettable via counter)

**Customization via `footnote.entry`:**
```typst
#show footnote.entry: set text(size: 0.9em)
```

**Example:**
```typst
Check the docs for details.#footnote[https://typst.app/docs]

#set footnote(numbering: "a")
Later footnotes use letter numbering.#footnote[Another note]
```

---

## Styling with Set & Show Rules

### set Rules

Set rules configure element properties persistently within their scope:

```typst
#set element(parameter: value)
```

**Characteristics:**
- Only modify optional parameters
- Scope: Top-level rules persist to end of file; nested rules scope to block
- Apply to all future instances of the element

**Examples:**
```typst
#set heading(numbering: "I.")          // All headings use Roman numerals
#set text(size: 12pt)                   // All text defaults to 12pt
#set par(justify: true)                 // All paragraphs justified
```

### show Rules

Show rules redefine element appearance via selectors and transformations.

**Form 1: Selector with set rule**
```typst
#show selector: set element(parameter: value)
```

**Form 2: Selector with transformation function**
```typst
#show selector: it => [custom formatting]
```

The parameter `it` references the matched element.

### Selector Patterns

- **Elements:** `show heading: ..` — all headings
- **Text strings:** `show "Text": ..` — literal text matches
- **Regex:** `show regex("\w+"): ..` — pattern-based
- **Conditional:** `show heading.where(level: 1): ..` — property filtering
- **Labels:** `show <label>: ..` — labeled content
- **Global:** `show: rest => ..` — all subsequent content

### Composition Example

```typst
#set heading(numbering: "1.")
#show heading: set align(center)
#show heading: it => block(
  above: 1em,
  below: 1em,
  text(size: 1.2em, weight: "bold")[#it.body]
)
```

Rules apply in order; later rules override earlier styling.

---

## Universe Packages

### ctheorems (v1.1.3)

**Purpose:** Numbered theorem environments with custom formatting and counter management.

**Import:**
```typst
#import "@preview/ctheorems:1.1.3": *
#show: thmrules
```

**Key Functions:**

- `thmbox(identifier, name, fill, inset, ...)` — Boxed theorem with background color
- `thmplain(identifier, name, ...)` — Plain (non-boxed) theorem
- `thmproof(name, ..., qed-symbol)` — Proof block with optional QED marker

**Core Parameters:**
- `identifier`: Counter name (shared across related environments)
- `base`: Attach counter to heading or element
- `base_level`: Manual depth specification
- `fill`: Background color
- `inset`: Internal padding
- `qed-symbol`: Proof ending marker (custom symbol)
- `numbering`: Counter display format (e.g., `"1.1"`)

**Example:**
```typst
#let theorem = thmbox(
  "theorem",
  "Theorem",
  fill: rgb("#eeffee"),
  inset: 1em,
)

#let lemma = thmbox(
  "lemma",
  "Lemma",
  fill: rgb("#ffffee"),
  base: "theorem",  // Share counter with theorem
)

#theorem[
  Statement here.
]

#lemma[
  Supporting statement.
]
```

**Cross-referencing:**
```typst
#let theorem = thmbox(...) <my-theorem>

See @my-theorem for details.
```

**License:** MIT

---

### scienceicons (v0.1.0)

**Purpose:** 29 SVG icons for open-science publishing (ORCID, GitHub, ArXiv, Creative Commons, etc.).

**Import and Usage:**
```typst
#import "@preview/scienceicons:0.1.0": *

This is open access #open-access-icon(height: 1.1em, baseline: 20%)
```

**Customization Parameters:**
- `color`: Color override (default: black)
- `height`: Icon size (default: `1.1em`)
- `baseline`: Vertical alignment (default: `13.5%`)

**Available Icons (29 total):**

**Academic/Research:**
- `arxiv-icon`, `orcid-icon`, `osi-icon`, `ror-icon`, `open-access-icon`

**Publishing Platforms:**
- `jupyter-icon`, `jupyter-book-icon`, `jupyter-text-icon`, `myst-icon`, `curvenote-icon`, `binder-icon`

**Creative Commons:**
- `cc-icon`, `cc-by-icon`, `cc-nc-icon`, `cc-nd-icon`, `cc-sa-icon`, `cc-zero-icon`

**Communication & Social:**
- `github-icon`, `discord-icon`, `discourse-icon`, `email-icon`, `linkedin-icon`, `mastodon-icon`, `slack-icon`, `twitter-icon`, `x-icon`, `bluesky-icon`, `youtube-icon`, `website-icon`

**Example:**
```typst
#let icon(name, color: black, height: 1em) = {
  let icon_fn = eval("@preview/scienceicons:0.1.0::" + name)
  icon_fn(color: color, height: height, baseline: 20%)
}

ORCID: #icon("orcid-icon", color: green, height: 1.2em)
GitHub: #icon("github-icon")
```

**License:** MIT | **Author:** rowanc1

---

## Quick Reference: Achieving Common Tasks

### Create a Two-Column Layout
```typst
#set page(columns: 2)
```

### Add Page Numbers
```typst
#set page(numbering: "1 / 1", number-align: center + bottom)
```

### Center Headings
```typst
#show heading: set align(center)
```

### Add First-Line Indent
```typst
#set par(first-line-indent: 1em)
```

### Justify Text
```typst
#set par(justify: true)
```

### Create Bold Headers
```typst
#set text(font: "Libertinus Serif")
#show heading: set text(weight: "bold", size: 14pt)
```

### Import and Use a Template
```typst
#import "template.typ": my_template
#show: my_template.with(title: "My Document", author: "Me")
```

### Customize Bibliography Style
```typst
#bibliography("refs.bib", style: "apa", full: false)
```

### Add Footnotes
```typst
Check this.#footnote[Source: example.com]
```

### Create Theorem Environments
```typst
#import "@preview/ctheorems:1.1.3": *
#show: thmrules
#let theorem = thmbox("thm", "Theorem", fill: rgb("#eee"))
#theorem[Statement]
```

### Display Science Icons
```typst
#import "@preview/scienceicons:0.1.0": orcid-icon, github-icon
ORCID: #orcid-icon(height: 1em) \
GitHub: #github-icon(height: 1em)
```

---

## Notes & Known Limitations

1. **PGF/TikZ Equivalent:** No direct equivalent yet; `cetz` package is the community alternative
2. **Binding Margins:** Use `margin: (inside: Xcm, outside: Ycm)` instead of `geometry`'s `twoside` margin differences
3. **Page Breaks:** Use `page()` function (causes explicit page break); negative padding with `pad()` for one-off margins
4. **Math Delimiters:** Auto-scale without explicit `\left`/`\right`
5. **Multi-char Subscripts:** Use parentheses: `$x_(a -> b)$`
6. **Floating Placement:** `placement: top` or `bottom` within `figure()`; scope via `scope: "column"` or `"parent"`

---

**End of Reference Document**


---
title: "Emoji Font Test"
author: "compiler diagnostic"
date: "2026-09-15"
tags: [test, emoji, compile-all]
---

**Version 1** (September 2026)

Raw emoji direct input, four codepoints (as used in the assessment status matrix):

✅ 🟠 🔴 🔥

If all four render in color below the heading, the TEXMFHOME emoji fix
(\texttt{compile-all.zsh} installing Noto Color Emoji and calling it by
family name) works. If they render as empty squares, the \texttt{emoji}
package still cannot use the font; in that case the PNG route (as in
\texttt{DOCS\_ASSESSMENT.md}) remains the fallback. Delete this file
once confirmed.

Run with:
`zsh /Users/joko/Git/ai/IT/scripts/compile-all.zsh "/Users/joko/Git/BFRep/(3)BeyondHulten/docs/archive/emoji_test.md"`

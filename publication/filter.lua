-- filter.lua — turn marked-up Divs into journal-quality LaTeX floats.
-- .figure / .figurewide -> figure / figure* ; .table -> table* with
-- proportional p{} columns ; .bibliography -> thebibliography (refNN keys).
--
-- pandoc.write only serialises the filter's OWN document, so Div handlers
-- stash closures and a final Pandoc pass renders them by swapping blocks.

local STASH = {}   -- list of function(doc) -> latex string

local function inlines_of(blocks)
  local out = {}
  for _, b in ipairs(blocks) do
    if b.t == 'Para' or b.t == 'Plain' then
      for _, i in ipairs(b.content) do table.insert(out, i) end
      table.insert(out, pandoc.SoftBreak())
    end
  end
  return out
end

local function all_rows(tbl)
  local rows = {}
  for _, row in ipairs(tbl.head.rows) do table.insert(rows, row) end
  for _, sec in ipairs(tbl.bodies) do
    for _, row in ipairs(sec.head) do table.insert(rows, row) end
    for _, row in ipairs(sec.body) do table.insert(rows, row) end
  end
  return rows
end

local function body_rows(tbl)
  local rows = {}
  for _, sec in ipairs(tbl.bodies) do
    for _, row in ipairs(sec.head) do table.insert(rows, row) end
    for _, row in ipairs(sec.body) do table.insert(rows, row) end
  end
  return rows
end

local RAWIFY = {
  Code = function(c)
    if c.text:match('^\\cite{') or c.text:match('^\\ref{')
       or c.text:match('^Figs?%.~\\ref{') or c.text:match('^Table~\\ref{')
       or c.text:match('^\\url{') then
      return pandoc.RawInline('latex', c.text)
    end
    return c
  end
}

local rawify

local function breakable_tt(s)
  -- escape LaTeX specials first, then allow breaks after / _ . separators
  s = s:gsub('\\', '\\textbackslash{}')
  s = s:gsub('([%&%%#_{}])', '\\%1')
  s = s:gsub('([/%_%.])', '%1\\allowbreak{}')
  return '\\texttt{' .. s .. '}'
end

local function rawify_one(i)
  if i.t == 'Code' and (i.text:match('^\\cite{') or i.text:match('^\\ref{')
      or i.text:match('^Figs?%.~\\ref{') or i.text:match('^Table~\\ref{') or i.text:match('^\\url{')) then
    return pandoc.RawInline('latex', i.text)
  end
  if i.t == 'Code' and #i.text > 28 and i.text:find('/') then
    return pandoc.RawInline('latex', breakable_tt(i.text))
  end
  for _, k in ipairs{'Strong', 'Emph', 'Superscript', 'Subscript', 'Span', 'Quoted'} do
    if i.t == k then
      local copy = i:clone()
      copy.content = rawify(copy.content)
      return copy
    end
  end
  if i.t == 'Link' then
    local copy = i:clone()
    copy.content = rawify(copy.content)
    return copy
  end
  return i
end

rawify = function(inlines)
  local out = {}
  for _, i in ipairs(inlines) do table.insert(out, rawify_one(i)) end
  return out
end

local function stash(fn)
  table.insert(STASH, fn)
  return string.format('\x01STASH%d\x01', #STASH)
end

-- DOI table: publication/doi.json (verified by fetch_dois.py against Crossref)
local DOIS = {}
do
  local f = io.open('doi.json', 'r')
  if f then
    local s = f:read('*a'); f:close()
    for k, v in s:gmatch('"(%d+)"%s*:%s*"([^"]+)"') do DOIS[tonumber(k)] = v end
  end
end

function Div(el)
  if el.classes:includes('figure') or el.classes:includes('figurewide') then
    local wide = el.classes:includes('figurewide')
    local id = el.identifier ~= '' and el.identifier or (el.attributes.id or '')
    local src = el.attributes.src or ''
    local blocks = el.content
    local tok = stash(function(doc)
      doc.blocks = pandoc.Blocks{ pandoc.Para(rawify(inlines_of(blocks))) }
      local cap = pandoc.write(doc, 'latex'):gsub('\n+$', '')
      return string.format(
        '\\begin{%s}[%s]\n\\centering\n\\includegraphics[width=%s]{%s}\n\\caption{%s}\n\\label{%s}\n\\end{%s}',
        wide and 'figure*' or 'figure',
        wide and 'ht' or 'ht',
        wide and '\\textwidth' or '\\columnwidth',
        src, cap, id, wide and 'figure*' or 'figure')
    end)
    return pandoc.RawBlock('latex', tok)
  end

  if el.classes:includes('table') then
    local id = el.identifier ~= '' and el.identifier or (el.attributes.id or '')
    local title, tbl = nil, nil
    for _, b in ipairs(el.content) do
      if b.t == 'Para' and not title then title = b.content end
      if b.t == 'Table' then tbl = b end
    end
    if not tbl then return el end
    local title_inlines = title or {}

    -- column widths proportional to the longest cell text
    local rows = all_rows(tbl)
    local ncol = #tbl.head.rows[1].cells
    local maxlen = {}
    for j = 1, ncol do maxlen[j] = #pandoc.utils.stringify(tbl.head.rows[1].cells[j]) end
    for _, row in ipairs(rows) do
      for j, cell in ipairs(row.cells) do
        maxlen[j] = math.max(maxlen[j], #pandoc.utils.stringify(cell))
      end
    end
    local total = 0
    for j = 1, ncol do total = total + maxlen[j] end
    -- longest single WORD per column: p{} columns cannot shrink below it
    local maxword = {}
    for j = 1, ncol do maxword[j] = 0 end
    for _, row in ipairs(rows) do
      for j, cell in ipairs(row.cells) do
        for word in (pandoc.utils.stringify(cell) .. ' '):gmatch('(%S+)%s') do
          maxword[j] = math.max(maxword[j], #word)
        end
      end
    end
    local frac = {}
    for j = 1, ncol do frac[j] = math.max(0.07, maxlen[j] / total) end
    local s = 0
    for j = 1, ncol do s = s + frac[j] end
    local f = {}
    for j = 1, ncol do f[j] = frac[j] / s * 0.97 end
    -- footnotesize ~4.1pt/char; +6pt cell padding; textwidth ~504pt
    local fmin = {}
    local sum, sumflex = 0, 0
    for j = 1, ncol do
      fmin[j] = (maxword[j] * 4.1 + 6) / 504
      f[j] = math.max(f[j], fmin[j])
      sum = sum + f[j]
    end
    -- if the floors pushed the total over budget, shave the flexible columns
    if sum > 0.97 then
      local excess = sum - 0.97
      for j = 1, ncol do sumflex = sumflex + math.max(0, f[j] - fmin[j]) end
      if sumflex > 0 then
        for j = 1, ncol do
          f[j] = math.max(fmin[j], f[j] - excess * math.max(0, f[j] - fmin[j]) / sumflex)
        end
      end
    end
    local cols = {}
    for j = 1, ncol do
      cols[j] = string.format('p{%.3f\\textwidth}', f[j])
    end

    local tok = stash(function(doc)
      local function latex_of(inlines)
        doc.blocks = pandoc.Blocks{ pandoc.Para(rawify(inlines)) }
        return pandoc.write(doc, 'latex'):gsub('\n+$', '')
      end
      local function rowlatex(cells, bold)
        local parts = {}
        for _, cell in ipairs(cells) do
          local t = latex_of(inlines_of(cell.contents))
          if bold then t = '\\textbf{' .. t .. '}' end
          -- p{} cells are justified by default: in narrow columns LaTeX
          -- stretches word spaces to the right margin ("étalé"). Force
          -- ragged-right inside each cell (array package is unavailable
          -- with revtex4-2, so wrap the content instead of the column).
          -- The explicit \par is required: the cell's implicit \par fires
          -- after the group closes, which would discard \raggedright.
          table.insert(parts, '{\\raggedright ' .. t .. '\\par}')
        end
        return table.concat(parts, ' & ') .. ' \\\\ \\addlinespace[2pt]'
      end
      local body = {}
      for _, row in ipairs(body_rows(tbl)) do table.insert(body, rowlatex(row.cells, false)) end
      return string.format(
        '\\begin{table*}[t]\n\\caption{%s}\n\\label{%s}\n\\footnotesize\n\\begin{tabular}{@{}%s@{}}\n\\toprule\n%s \\midrule\n%s \\bottomrule\n\\end{tabular}\n\\end{table*}',
        latex_of(title_inlines), id, table.concat(cols, ''),
        rowlatex(tbl.head.rows[1].cells, true),
        table.concat(body, '\n\\grayhline\n'))
    end)
    return pandoc.RawBlock('latex', tok)
  end

  if el.classes:includes('bibliography') then
    local ol = el.content[1]
    if not ol or ol.t ~= 'OrderedList' then return el end
    local tok = stash(function(doc)
      local items = {}
      for i, blocks in ipairs(ol.content) do
        doc.blocks = pandoc.Blocks{ pandoc.Para(rawify(inlines_of(blocks))) }
        local t = pandoc.write(doc, 'latex'):gsub('\n+$', '')
        if DOIS[i] then
          t = t .. ' \\url{https://doi.org/' .. DOIS[i] .. '}'
        end
        table.insert(items, string.format('\\bibitem{ref%02d} %s', i, t))
      end
      return '\\begin{thebibliography}{99}\n' .. table.concat(items, '\n\n') ..
             '\n\\end{thebibliography}'
    end)
    return pandoc.RawBlock('latex', tok)
  end

  return el
end

function Pandoc(doc)
  -- backtick spans holding LaTeX (citations, cross-refs) -> RawInline
  doc = doc:walk({
    Code = function(c)
      if c.text:match('^\\cite{') or c.text:match('^\\ref{')
         or c.text:match('^Figs?%.~\\ref{') or c.text:match('^Table~\\ref{')
         or c.text:match('^\\url{') then
        return pandoc.RawInline('latex', c.text)
      end
      if #c.text > 28 and c.text:find('/') then
        return pandoc.RawInline('latex', breakable_tt(c.text))
      end
      return c
    end
  })
  if #STASH == 0 then return doc end
  local rendered = {}
  for i, fn in ipairs(STASH) do
    local saved = doc.blocks
    rendered[i] = fn(doc)
    doc.blocks = saved
  end
  doc = doc:walk({
    RawBlock = function(rb)
      if rb.format ~= 'latex' then return rb end
      local out = rb.text:gsub('\x01STASH(%d+)\x01', function(n)
        return rendered[tonumber(n)] or ''
      end)
      return pandoc.RawBlock('latex', out)
    end
  })
  return doc
end

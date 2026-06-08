if not FORMAT:match("latex") then
  return
end

function Table(table)
  if table.classes:includes("format-table") then
    return {
      pandoc.RawBlock("latex", "\\begingroup\\small"),
      table,
      pandoc.RawBlock("latex", "\\endgroup"),
    }
  end
end

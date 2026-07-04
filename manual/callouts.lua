-- Map fenced-div classes to LaTeX callout environments (manual PDF only).
--   ::: note      -> notebox
--   ::: warning   -> warnbox
--   ::: tip       -> tipbox
local map = { note = "notebox", warning = "warnbox", tip = "tipbox" }

function Div(el)
  for _, cls in ipairs(el.classes) do
    local env = map[cls]
    if env then
      local out = pandoc.List()
      out:insert(pandoc.RawBlock("latex", "\\begin{" .. env .. "}"))
      out:extend(el.content)
      out:insert(pandoc.RawBlock("latex", "\\end{" .. env .. "}"))
      return out
    end
  end
  return nil
end

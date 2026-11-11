/**
 * Compact contextual hint with an accessible hover/focus tooltip.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { HelpCircle } from "lucide-react";

export function InfoHint({ text }: { text: string }) {
  return (
    <span className="info-hint" tabIndex={0} aria-label={text}>
      <HelpCircle size={13} />
      <span className="info-hint-tooltip" role="tooltip">
        {text}
      </span>
    </span>
  );
}

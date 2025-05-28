// docs/source/_static/custom.js
document.addEventListener('DOMContentLoaded', function() {
    // Nuclear option - remove all commas and clean whitespace
    const signatures = document.querySelectorAll('.sig');
    signatures.forEach(sig => {
        // Remove all comma elements
        const commas = sig.querySelectorAll('.sig-paren');
        commas.forEach(comma => comma.remove());
        
        // Remove any text nodes containing commas
        const walker = document.createTreeWalker(sig, NodeFilter.SHOW_TEXT);
        const nodes = [];
        while (walker.nextNode()) {
            if (walker.currentNode.nodeValue.includes(',')) {
                nodes.push(walker.currentNode);
            }
        }
        nodes.forEach(node => node.remove());
        
        // Clean up parameter spacing
        const params = sig.querySelectorAll('.sig-param');
        params.forEach(param => {
            param.innerHTML = param.innerHTML.replace(/\s*,\s*/g, '');
        });
    });
});

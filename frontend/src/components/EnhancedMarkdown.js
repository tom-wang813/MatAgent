import React from 'react';
import ReactMarkdown from 'react-markdown';
import remarkMath from 'remark-math';
import remarkGfm from 'remark-gfm';
import rehypeKatex from 'rehype-katex';
import rehypeHighlight from 'rehype-highlight';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { vscDarkPlus } from 'react-syntax-highlighter/dist/esm/styles/prism';
import { Box, Typography } from '@mui/material';
import 'katex/dist/katex.min.css';
import 'highlight.js/styles/github.css'; // GitHub style for code highlighting

const EnhancedMarkdown = ({ children, ...props }) => {
  const components = {
    // Code blocks with syntax highlighting
    code({ node, inline, className, children, ...props }) {
      const match = /language-(\w+)/.exec(className || '');
      const language = match ? match[1] : '';
      
      if (!inline && language) {
        return (
          <Box sx={{ my: 2 }}>
            <SyntaxHighlighter
              style={vscDarkPlus}
              language={language}
              PreTag="div"
              {...props}
            >
              {String(children).replace(/\n$/, '')}
            </SyntaxHighlighter>
          </Box>
        );
      }
      
      // Inline code
      return (
        <Box
          component="code"
          sx={{
            backgroundColor: 'grey.100',
            color: 'error.main',
            px: 0.5,
            py: 0.25,
            borderRadius: 0.5,
            fontFamily: 'monospace',
            fontSize: '0.875em',
          }}
          {...props}
        >
          {children}
        </Box>
      );
    },
    
    // Enhanced headings
    h1: ({ children }) => (
      <Typography variant="h4" component="h1" sx={{ mt: 3, mb: 2, fontWeight: 'bold' }}>
        {children}
      </Typography>
    ),
    h2: ({ children }) => (
      <Typography variant="h5" component="h2" sx={{ mt: 2.5, mb: 1.5, fontWeight: 'bold' }}>
        {children}
      </Typography>
    ),
    h3: ({ children }) => (
      <Typography variant="h6" component="h3" sx={{ mt: 2, mb: 1, fontWeight: 'bold' }}>
        {children}
      </Typography>
    ),
    
    // Enhanced paragraphs
    p: ({ children }) => (
      <Typography variant="body1" component="p" sx={{ mb: 1.5, lineHeight: 1.6 }}>
        {children}
      </Typography>
    ),
    
    // Enhanced lists
    ul: ({ children }) => (
      <Box component="ul" sx={{ pl: 2, mb: 1.5 }}>
        {children}
      </Box>
    ),
    ol: ({ children }) => (
      <Box component="ol" sx={{ pl: 2, mb: 1.5 }}>
        {children}
      </Box>
    ),
    li: ({ children }) => (
      <Typography component="li" variant="body1" sx={{ mb: 0.5 }}>
        {children}
      </Typography>
    ),
    
    // Enhanced blockquotes
    blockquote: ({ children }) => (
      <Box
        sx={{
          borderLeft: '4px solid',
          borderColor: 'primary.main',
          backgroundColor: 'grey.50',
          pl: 2,
          py: 1,
          my: 2,
          fontStyle: 'italic',
        }}
      >
        {children}
      </Box>
    ),
    
    // Enhanced tables
    table: ({ children }) => (
      <Box
        sx={{
          width: '100%',
          overflowX: 'auto',
          my: 2,
          border: '1px solid',
          borderColor: 'divider',
          borderRadius: 1,
        }}
      >
        <Box component="table" sx={{ width: '100%', borderCollapse: 'collapse' }}>
          {children}
        </Box>
      </Box>
    ),
    thead: ({ children }) => (
      <Box component="thead" sx={{ backgroundColor: 'grey.100' }}>
        {children}
      </Box>
    ),
    th: ({ children }) => (
      <Typography
        component="th"
        variant="subtitle2"
        sx={{
          p: 1.5,
          textAlign: 'left',
          borderBottom: '1px solid',
          borderColor: 'divider',
          fontWeight: 'bold',
        }}
      >
        {children}
      </Typography>
    ),
    td: ({ children }) => (
      <Typography
        component="td"
        variant="body2"
        sx={{
          p: 1.5,
          borderBottom: '1px solid',
          borderColor: 'divider',
        }}
      >
        {children}
      </Typography>
    ),
    
    // Enhanced links
    a: ({ href, children }) => (
      <Box
        component="a"
        href={href}
        target="_blank"
        rel="noopener noreferrer"
        sx={{
          color: 'primary.main',
          textDecoration: 'none',
          '&:hover': {
            textDecoration: 'underline',
          },
        }}
      >
        {children}
      </Box>
    ),
    
    // Enhanced horizontal rules
    hr: () => (
      <Box
        sx={{
          border: 'none',
          height: '1px',
          backgroundColor: 'divider',
          my: 3,
        }}
      />
    ),
  };

  return (
    <ReactMarkdown
      remarkPlugins={[remarkMath, remarkGfm]}
      rehypePlugins={[rehypeKatex, rehypeHighlight]}
      components={components}
      {...props}
    >
      {children}
    </ReactMarkdown>
  );
};

export default EnhancedMarkdown;

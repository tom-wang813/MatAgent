import React, { useState } from 'react';
import { Box, Typography, Accordion, AccordionSummary, AccordionDetails, Button, CircularProgress, Link, Paper } from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import DownloadIcon from '@mui/icons-material/Download';
import SettingsIcon from '@mui/icons-material/Settings';
import CheckCircleIcon from '@mui/icons-material/CheckCircle';
import CancelIcon from '@mui/icons-material/Cancel';
import AccessTimeIcon from '@mui/icons-material/AccessTime';

const ToolCallDisplay = ({ toolCalls = [] }) => {
  const [expanded, setExpanded] = useState({});

  const handleChange = (panel) => (event, isExpanded) => {
    setExpanded((prev) => ({
      ...prev,
      [panel]: isExpanded ? true : false,
    }));
  };

  const renderToolOutput = (output) => {
    try {
      const parsedOutput = JSON.parse(output);
      
      if (parsedOutput.file_path) {
        const downloadUrl = `/api/download/mcp_file?file_path=${encodeURIComponent(parsedOutput.file_path)}`;
        return (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mt: 1 }}>
            {parsedOutput.message && <Typography variant="body2">{parsedOutput.message}</Typography>}
            <Button
              variant="outlined"
              startIcon={<DownloadIcon />}
              href={downloadUrl}
              target="_blank"
              rel="noopener noreferrer"
              size="small"
            >
              下载文件: {parsedOutput.file_path.split('/').pop()}
            </Button>
          </Box>
        );
      } else if (parsedOutput.svg_data) {
        return (
          <Box
            sx={{
              mt: 1,
              '& svg': {
                maxWidth: '100%',
                height: 'auto',
                display: 'block',
                margin: 'auto',
              },
            }}
            dangerouslySetInnerHTML={{ __html: parsedOutput.svg_data }}
          />
        );
      } else if (parsedOutput.png_base64) {
        return (
          <Box sx={{ mt: 1, display: 'flex', justifyContent: 'center' }}>
            <img src={`data:image/png;base64,${parsedOutput.png_base64}`} alt="Molecule" style={{ maxWidth: '100%', height: 'auto' }} />
          </Box>
        );
      } else {
        return (
          <Typography component="pre" sx={{ bgcolor: 'grey.100', p: 1, borderRadius: 1, overflowX: 'auto', mt: 1 }}>
            {JSON.stringify(parsedOutput, null, 2)}
          </Typography>
        );
      }
    } catch (e) {
      return (
        <Typography component="pre" sx={{ bgcolor: 'grey.100', p: 1, borderRadius: 1, overflowX: 'auto', mt: 1 }}>
          {output}
        </Typography>
      );
    }
  };

  if (!toolCalls || toolCalls.length === 0) {
    return null;
  }

  return (
    <Box sx={{ mt: 2, mb: 2 }}>
      <Typography variant="h6" sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 1 }}>
        <SettingsIcon fontSize="small" />
        工具调用 ({toolCalls.length})
      </Typography>
      
      {toolCalls.map((call, index) => (
        <Paper key={index} elevation={1} sx={{ mb: 1 }}>
          <Accordion expanded={expanded[`panel${index}`]} onChange={handleChange(`panel${index}`)}>
            <AccordionSummary
              expandIcon={<ExpandMoreIcon />}
              aria-controls={`panel${index}bh-content`}
              id={`panel${index}bh-header`}
              sx={{ bgcolor: 'grey.50' }}
            >
              <Box sx={{ display: 'flex', alignItems: 'center', width: '100%', pr: 2 }}>
                <Typography variant="subtitle1" sx={{ flexGrow: 1 }}>
                  {call.function.name}
                </Typography>
                {call.status === 'success' && <CheckCircleIcon color="success" fontSize="small" />}
                {call.status === 'error' && <CancelIcon color="error" fontSize="small" />}
                {call.status === 'running' && <CircularProgress size={16} sx={{ color: 'text.secondary' }} />}
                {call.status === 'pending' && <AccessTimeIcon color="action" fontSize="small" />}
              </Box>
            </AccordionSummary>
            <AccordionDetails>
              <Box sx={{ mb: 2 }}>
                <Typography variant="subtitle2" gutterBottom>输入参数:</Typography>
                <Typography component="pre" sx={{ bgcolor: 'grey.100', p: 1, borderRadius: 1, overflowX: 'auto' }}>
                  {JSON.stringify(JSON.parse(call.function.arguments), null, 2)}
                </Typography>
              </Box>
              {call.output && (
                <Box>
                  <Typography variant="subtitle2" gutterBottom>输出结果:</Typography>
                  {renderToolOutput(call.output)}
                </Box>
              )}
            </AccordionDetails>
          </Accordion>
        </Paper>
      ))}
    </Box>
  );
};

export default ToolCallDisplay;

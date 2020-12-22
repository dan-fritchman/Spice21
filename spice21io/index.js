
// 
// # Spice21 Simulation Server 
// 

const spice21 = require('../spice21js');

const bodyParser = require('body-parser');
const express = require('express');

const CONTENT_TYPE = 'application/protobuf';

const app = express();
app.use(bodyParser.raw({ type: CONTENT_TYPE, limit: '2mb' }))

app.get('/', (_req, res) => { 
  res.send(`
    Thanks for checking out Spice21. 
    (You probably want to POST here rather than GET!)
    For the latest docs check out https://github.com/HW21/Spice21
  `);
});

app.post('/echo', (req, res) => {
  // Echo the request back 
  res.setHeader('content-type', CONTENT_TYPE);
  res.send(req.body);
});

const port = process.env.PORT || 8080;
app.listen(port, () => {
  console.log(`Server up, listening on port ${port}`);
});

// Exports for testing purposes.
module.exports = app;

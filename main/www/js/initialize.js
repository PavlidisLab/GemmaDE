var fileValidated = false;

$(document).one('shiny:connected', function() {
  Shiny.addCustomMessageHandler('querySet', querySet);
  Shiny.addCustomMessageHandler('queryReset', queryReset);
  Shiny.addCustomMessageHandler('fileUpload', message => fileValidated = message);
  Shiny.setInputValue('LOAD', true);
  
  // Clickable examples below entries
  loadExamples();
});

$(document).on('shiny:value', function(event) {
  if(event.name === 'results_genes') {
    setTimeout(function() {
      updatePanels();
    }, 100);
  } else if(event.name == 'results_header') {
    Shiny.setInputValue('UPDATED', Math.random());
  }
});

$(document).on('shiny:inputchanged', function(event) {
  if(event.name === 'genes')
    validate();
});

// Make the examples clickable
function loadExamples() {
  $('a[genes]').click(function() {
    var sz = $('#genes').selectize()[0].selectize;
    
    sz.clear(true);
    
    if(JSON.parse($(this).attr('genes')) == 'random()') {
      Shiny.setInputValue('RANDOM_GENES', Math.random());
    } else {
      addGenes(JSON.parse($(this).attr('genes')));
    }
  });
}

function queryReset(genes) {
  querySet(genes, true);
}

// Forcibly input genes from a list.
function querySet(genes, clear = false) { // TODO This is hacky. Might not work
  setTimeout(function() { addGenes(genes, clear); }, 100);
}

// Add genes to the genes bar
function addGenes(genes, clear = false) {
  var sz = $('#genes').selectize()[0].selectize;
  
  if(clear)
    sz.clear(true);
  
  if(genes instanceof Array)
    genes.forEach(gene => sz.createItem(gene, false));
  else
    sz.createItem(genes, false);
      
  Shiny.setInputValue('genes', genes);
}

// Disable the search button if the gene query is empty and no file is uploaded
function validate() {
  var Ng = $('#genes').selectize()[0].selectize.getValue().length;
  var accept = Ng !== 0 || fileValidated;
  
  $('#search').prop('disabled', !accept);
}

function onTableCreated() {
  $('[data-toggle="tooltip"]').tooltip();
}

$(document).click(function (e) {
  if($(e.target).parent().find('[data-toggle="popover"]').length > 0) {
    $('[data-toggle="popover"]').popover('hide');
    e.preventDefault();
  }
});

function addMathJax() {
  if(window.MathJax) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
}

function onTableDraw() {
  $('[data-toggle="popover"]').popover();
  
  $('[data-toggle="popover"]').click(function(e) {
    e.preventDefault();
    
    $('[data-toggle="popover"]').not(this).popover('hide');
    $(this).popover('toggle');
  });
  
  addMathJax();
  
  $('.spark:not(:has(canvas))').each(function(index) {
      $(this).sparkline('html', {
        type: $(this).attr('type'),
        sliceColors: ['#4DAC26', '#D7191C'],
        barColor: '#002145',
        zeroColor: '#002145',
        chartRangeMin: 0,
        chartRangeMax: 1
      });
  });
}

function histogram(data) {
  var bins = [];
  
  for(i = 0.1; i <= 1; i += 0.1) {
    bins.push(data.filter(x => x >= (i - 0.1) && x < i).length);
  }
  
  return bins;
}

const mean = arr => arr.reduce((sume, el) => sume + el, 0) / arr.length;

function asSparkline2(data, type, row, meta) {
  data = data.split(',').map(x => +x);
  
  return type === 'display' ?
    '<span class="spark" type="pie" style="display: inline-block; width: 100%;" mean="' +
      mean(data) + '">' + data + '</span>' :
    data;
}

function asSparkline(data, type, row, meta) {
  data = data.split(',').map(x => +x);
  
  return type === 'display' ?
    '<span class="spark" type="bar" style="display: inline-block; width: 100%;" mean="' +
      mean(data) + '">' + histogram(data) + '</span>' :
    data;
}

function asPval(data, type, row, meta) {
  if(data === null) return data;
  
  if(Math.floor(Math.log10(data)) >= -2) pval = Number(Math.round(data + 'e3') + 'e-3');
  else pval = data.toExponential(3).replace('e', ' * 10^{') + '}';
  
  return type === 'display' ?
    '<span class="pvalue" pval="' + data + '">$$' + pval + '$$</span>' :
  data;
}

function unformatSpark(html, row, column, node) {
  if($(node).children().hasClass('spark'))
    return $(html).attr('mean');
  else if($(node).children().hasClass('pvalue'))
    return $(html).attr('pval');
  return html;
}

function plotData(e, dt, node, config) {
  Shiny.setInputValue('plotData', {
    page: dt.rows().indexes().toArray(),
    selected: dt.rows('.selected').indexes().toArray()
  }, { priority: 'event' });
}
